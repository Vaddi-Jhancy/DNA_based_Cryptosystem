import argparse
import os
from cryptography.hazmat.primitives import hashes, serialization
from cryptography.hazmat.primitives.asymmetric import ec, rsa, padding as rsa_padding
from cryptography.exceptions import InvalidSignature
from Crypto.Cipher import AES
from Crypto.Util.Padding import unpad
from reedsolo import RSCodec
import heapq
from collections import defaultdict
import time
import hashlib
import json

with open('dna_tables.json', 'r') as f1:
    dna_tables = json.load(f1)

# DNA tables
DNA_ENCODING_TABLE = dna_tables["DNA_ENCODING_TABLE"]
DELIMITER = dna_tables["DELIMITER"]
# 64 codons for 64 possible values (6 bits per codon)
CODONS = dna_tables["CODONS"]

DNA_DECODING_TABLE = {v: k for k, v in DNA_ENCODING_TABLE.items()}

# Table 2 as a dictionary
table2_substitution = {
    tuple(k.split(',')): v
    for k, v in dna_tables["TABLE2_SUBSTITUTION"].items()
}

reverse_table2_substitution = {}

for (msg_char, cover_char), result_char in table2_substitution.items():
    # Key is (result_char, cover_char)
    reverse_table2_substitution[(result_char, cover_char)] = msg_char


def binary_to_codons(binary_str):
    # Pad to multiple of 6 bits
    if len(binary_str) % 6 != 0:
        binary_str += '0' * (6 - len(binary_str) % 6)
    codons = []
    for i in range(0, len(binary_str), 6):
        idx = int(binary_str[i:i+6], 2)
        codons.append(CODONS[idx])
    return ''.join(codons)

def codons_to_binary(codon_str):
    codon_map = {codon: format(i, '06b') for i, codon in enumerate(CODONS)}
    binary = ''
    for i in range(0, len(codon_str), 3):
        codon = codon_str[i:i+3]
        binary += codon_map[codon]
    return binary

def find_adjacent_repeats_indices(dna_sequence):
    indices = []
    for i in range(1, len(dna_sequence)):
        if dna_sequence[i] == dna_sequence[i - 1]:
            indices.append(i)
    return indices



# File operations
def read_file(path, binary=False):
    with open(path, 'rb' if binary else 'r') as f:
        return f.read().strip()

# DNA to binary
def dna_to_binary(dna_seq):
    return ''.join(DNA_DECODING_TABLE[nuc] for nuc in dna_seq)

def binary_to_dna(binary_sequence):
    return ''.join(DNA_ENCODING_TABLE[binary_sequence[i:i+2]] for i in range(0, len(binary_sequence), 2))

def xor_dna_sequences(dna1, dna2):
    # Convert DNA to binary
    bin1 = dna_to_binary(dna1)
    bin2 = dna_to_binary(dna2)
    # XOR the binaries
    xor_bin = ''.join(str(int(b1) ^ int(b2)) for b1, b2 in zip(bin1, bin2))
    # Convert back to DNA
    return binary_to_dna(xor_bin)

# Extract delimited binary
def extract_delimited_binary(dna_binary):
    start = dna_binary.find(DELIMITER)
    end = dna_binary.find(DELIMITER, start + len(DELIMITER))
    if start == -1 or end == -1 or start == end:
        raise ValueError("Delimiter not found properly.")
    return dna_binary[start + len(DELIMITER):end]

# Reed-Solomon decoding
def reed_solomon_decode(binary_str, ecc_symbols=10):
    byte_array = int(binary_str, 2).to_bytes((len(binary_str) + 7) // 8, 'big')
    rsc = RSCodec(ecc_symbols)
    # decoded = rsc.decode(byte_array)
    decoded = rsc.decode(byte_array)[0]
    return decoded

# Hamming decoding (12 bits to 8 bits)
def hamming_decode(data):
    decoded = ""
    for i in range(0, len(data), 12):
        block = data[i:i+12]
        if len(block) < 12:
            break
        d = [int(b) for b in block]
        # extract data bits
        data_bits = [d[2], d[4], d[5], d[6], d[8], d[9], d[10], d[11]]
        decoded += ''.join(str(b) for b in data_bits)
    return decoded

def parse_key(key_input):
    # If it's a hex string (64 characters), parse it
    if len(key_input) == 64 and all(c in "0123456789abcdefABCDEF" for c in key_input):
        return bytes.fromhex(key_input)

    # If it's a file path, try reading it as binary
    try:
        with open(key_input, "rb") as f:
            key = f.read()
            if len(key) != 32:
                raise ValueError(f"Key file must be 32 bytes. Found {len(key)} bytes.")
            return key
    except FileNotFoundError:
        raise ValueError("Key file not found.")
    except Exception as e:
        raise ValueError(f"Failed to load key: {e}")


# AES decryption (CBC)
def aes_decrypt_cbc(encrypted_data, key):
    iv = encrypted_data[:16]
    cipher = AES.new(key, AES.MODE_CBC, iv)
    decrypted = cipher.decrypt(encrypted_data[16:])
    return unpad(decrypted, AES.block_size)

# AES decryption (GCM)
def aes_decrypt_gcm(encrypted_data, key):
    nonce = encrypted_data[:16]
    tag = encrypted_data[16:32]
    ciphertext = encrypted_data[32:]
    cipher = AES.new(key, AES.MODE_GCM, nonce)
    try:
        decrypted_data = cipher.decrypt_and_verify(ciphertext, tag)
        print("✅ Integrity verified. Data has not been tampered with.")
        return decrypted_data
    except ValueError:
        print("❌ Integrity check failed! Data has been altered or is corrupt.")
        return None

# Huffman decompression
class HuffmanNode:
    def __init__(self, char=None, freq=0):
        self.char = char
        self.freq = freq
        self.left = None
        self.right = None

    def __lt__(self, other):
        return self.freq < other.freq

def rebuild_huffman_tree(codebook):
    root = HuffmanNode()
    for char, code in codebook.items():
        node = root
        for bit in code:
            if bit == '0':
                if not node.left:
                    node.left = HuffmanNode()
                node = node.left
            else:
                if not node.right:
                    node.right = HuffmanNode()
                node = node.right
        node.char = char
    return root

def huffman_decompress(binary_data, codebook, compressed_length):
    root = rebuild_huffman_tree(codebook)
    result = []
    node = root
    bit_count = 0
    for bit in binary_data:
        if bit_count >= compressed_length:
            break  # Stop once we reach the original compressed length
        node = node.left if bit == '0' else node.right
        if node.char is not None:
            result.append(node.char)
            node = root
        bit_count += 1
    return ''.join(result)

# Signature verification
def verify_signature(signature_scheme, public_key_data, signature_hex, hashed_data):
    signature = bytes.fromhex(signature_hex)

    if signature_scheme == 'rsa':
        public_key = serialization.load_ssh_public_key(public_key_data.encode())
        public_key.verify(
            signature,
            hashed_data,
            rsa_padding.PSS(mgf=rsa_padding.MGF1(hashes.SHA256()), salt_length=rsa_padding.PSS.MAX_LENGTH),
            hashes.SHA256()
        )
    else:
        public_key = serialization.load_pem_public_key(public_key_data.encode())
        public_key.verify(
            signature,
            hashed_data,
            ec.ECDSA(hashes.SHA256())
        )


def extract_delimited_encrypted_only(modified_binary_str):
    extracted_encrypted = []
    for i in range(2, len(modified_binary_str), 3):
        extracted_encrypted.append(modified_binary_str[i])
    return ''.join(extracted_encrypted)

# def extract_encrypted_dna_from_host(dna_sequence, repeated_nucleo_indexes):
#     return ''.join(dna_sequence[i] for i in repeated_nucleo_indexes)

def extract_encrypted_message_from_dna(encoded_dna, original_dna, repeated_indices):
    """
    Extracts the encrypted message embedded in the DNA using reverse of Table 2.
    """
    extracted = []
    for index in repeated_indices:
        embedded_char = encoded_dna[index]
        original_char = original_dna[index]
        # Lookup the reverse substitution
        msg_char = reverse_table2_substitution.get((embedded_char, original_char))
        if msg_char is not None:
            extracted.append(msg_char)
        else:
            # Something's wrong, fallback or handle silently
            pass
    return ''.join(extracted)

def hash_signature(signature_bytes):
    return hashlib.sha256(signature_bytes).hexdigest()

def dna_to_hex(dna_seq):
    """Convert DNA sequence back to hexadecimal string"""
    hex_str = []
    for i in range(0, len(dna_seq), 2):
        # Get two DNA nucleotides (4 bits total)
        pair = dna_seq[i:i+2]
        # Convert to 4-bit binary
        binary = DNA_DECODING_TABLE[pair[0]] + DNA_DECODING_TABLE[pair[1]]
        # Convert to hex character
        hex_str.append(format(int(binary, 2), '01x'))
    return ''.join(hex_str)

def main():

    with open('decode_config.json', 'r') as f:
        config = json.load(f)

    start_total = time.time()
    encoded_dna = read_file(config["encoded_dna"])
    dna_sequence = read_file(config["dna_sequence"])
    method = config["method"]
    dna_seq_2 = read_file(config["dna_file_2"]) if method == "method2" else None
    signature_hex = read_file(config["signature_file"])
    public_key = read_file(config["public_key_file"])
    with open(config["codebook_file"], "r") as f:
        codebook = json.load(f)
    aes_key = parse_key(config["keyfile"])
    cipher = config["cipher"]
    error_method = config["error"]
    sig_scheme = config.get("sig")
    compressed_length = int(read_file(config["compressed_length"]))

    # Extract and verify signature hash
    embedded_dna_hash = encoded_dna[10:30]  # Assuming 10-character hash embedded at position 10
    cleaned_encoded_dna = encoded_dna[:10] + encoded_dna[30:]  # Remove embedded hash for further processing
    signature_hash = dna_to_hex(embedded_dna_hash)  # Convert back to hex

    # Hash encoded DNA (without embedded hash)
    digest = hashes.Hash(hashes.SHA256())
    digest.update(cleaned_encoded_dna.encode())
    hashed = digest.finalize()

    # Step 1: Verify signature
    if sig_scheme is not None:
        start = time.time()
        try:
            verify_signature(sig_scheme, public_key, signature_hex, hashed)
            print("[✔] Signature verified.")
        except InvalidSignature:
            print("[✘] Signature verification failed.")
            return
        print(f"Signature Verification Time: {time.time() - start:.4f} seconds")

    # Verify that the signature hash was embedded correctly
    expected_hash = hash_signature(bytes.fromhex(signature_hex))[:10]
    if signature_hash != expected_hash:
        print("[✘] Embedded signature hash mismatch – potential tampering detected.")
        return
    print("[✔] Embedded signature hash verified.")


    start = time.time()
    if method=='method1':
    # Step 2: Decode and extract
        encoded_dna_binary = dna_to_binary(cleaned_encoded_dna)
        delimited_binary = extract_delimited_encrypted_only(encoded_dna_binary)
        encrypted_binary = extract_delimited_binary(delimited_binary)
        # print("encrypted_binary",encrypted_binary)

    else:
            # 1. Extract encrypted DNA from host sequence
        repeated_nucleo_indexes = find_adjacent_repeats_indices(dna_sequence)
        # print(f"repeated_indexes: {repeated_nucleo_indexes}")
        encrypted_dna = extract_encrypted_message_from_dna(cleaned_encoded_dna, dna_sequence, repeated_nucleo_indexes)
        
        # 2. Decrypt XOR
        decrypted_dna = xor_dna_sequences(encrypted_dna, dna_seq_2[:len(encrypted_dna)])
        # print("decrypted_dna",decrypted_dna)
        # print(f"Encrypted DNA length: {len(encrypted_dna)}")
        # print(f"Decrypted codon DNA length: {len(decrypted_dna)} (should be multiple of 3)")
        # print(f"Decrypted DNA (preview): {decrypted_dna[:30]}")

        
        # 3. Convert codons to binary
        delimited_binary = codons_to_binary(decrypted_dna)
        
        # 4. Remove delimiters
        encrypted_binary = extract_delimited_binary(delimited_binary)
    print(f"Encrypted Message Extraction Time: {time.time() - start:.4f} seconds")

    # Step 3: Error correction
    start = time.time()
    if error_method == "reed-solomon":
        encrypted_bytes = reed_solomon_decode(encrypted_binary)
    else:
        corrected_binary = hamming_decode(encrypted_binary)
        encrypted_bytes = int(corrected_binary, 2).to_bytes((len(corrected_binary) + 7) // 8, 'big')
    print(f"Error Correction Time: {time.time() - start:.4f} seconds")

    # Step 4: AES decryption
    start = time.time()
    if cipher == "cbc":
        decrypted_data = aes_decrypt_cbc(encrypted_bytes, aes_key)
    else:
        decrypted_data = aes_decrypt_gcm(encrypted_bytes, aes_key)
    print(f"AES Decryption Time: {time.time() - start:.4f} seconds")

    # Step 5: Huffman decompression
    start = time.time()
    binary_str = ''.join(format(byte, '08b') for byte in decrypted_data)
    # _>checked print(binary_str)
    # Truncate to the exact compressed length
    binary_str = binary_str[-compressed_length:] if len(binary_str) > compressed_length else binary_str
    original_message = huffman_decompress(binary_str, codebook, compressed_length)
    print(f"Huffman Decompression Time: {time.time() - start:.4f} seconds")

    print("\n[✓] Original Secret Message:")
    print(original_message)
    # print("compressed binary:",encrypted_binary)
    total_time = time.time() - start_total
    print(f"\nTotal Decryption Time: {total_time:.4f} seconds")

if __name__ == "__main__":
    main()
