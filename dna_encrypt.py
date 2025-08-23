import argparse
import time
import os
from cryptography.hazmat.primitives import serialization, hashes
from cryptography.hazmat.primitives.asymmetric import ec, rsa, padding
from Crypto.Cipher import AES
from Crypto.Util.Padding import pad, unpad
import heapq
from collections import defaultdict
from reedsolo import RSCodec
import json
import hashlib
import json

with open('dna_tables.json', 'r') as f1:
    dna_tables = json.load(f1)

DNA_ENCODING_TABLE = dna_tables["DNA_ENCODING_TABLE"]
# DELIMITER = dna_tables["DELIMITER"]
# 64 codons for 64 possible values (6 bits per codon)
CODONS = dna_tables["CODONS"]
DNA_DECODING_TABLE = {v: k for k, v in DNA_ENCODING_TABLE.items()}

# Convert string keys back to tuple keys for table2_substitution
table2_substitution = {
    tuple(k.split(',')): v
    for k, v in dna_tables["TABLE2_SUBSTITUTION"].items()
}

# Utility functions for file operations
def read_file(file_path):
    with open(file_path, 'r') as file:
        return file.read().strip()

def save_file(filename, data):
    with open(filename, "w") as file:
        file.write(data)
    print(f"Data successfully saved to {filename}")

def save_dic(file_path, data):
    with open(file_path,"w") as file:
        json.dump(data,file)

def add_delimiters(encrypted_binary, delimiter="10101010101010101010101010101010101010101010101010111111"):
    # encrypted_binary_str = ''.join(format(byte, '08b') for byte in encrypted_binary)
    # print(f"encrypted_binary: {encrypted_binary} ")
    return delimiter + encrypted_binary + delimiter

def dna_to_binary(dna_sequence):
    return ''.join(DNA_DECODING_TABLE[nucleotide] for nucleotide in dna_sequence)

def binary_to_dna(binary_sequence):
    return ''.join(DNA_ENCODING_TABLE[binary_sequence[i:i+2]] for i in range(0, len(binary_sequence), 2))

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

def xor_dna_sequences(dna1, dna2):
    # Convert DNA to binary
    bin1 = dna_to_binary(dna1)
    bin2 = dna_to_binary(dna2)
    # XOR the binaries
    xor_bin = ''.join(str(int(b1) ^ int(b2)) for b1, b2 in zip(bin1, bin2))
    # Convert back to DNA
    return binary_to_dna(xor_bin)

# def embed_in_second_dna(encrypted_dna, cover_dna):
#     # Replace the second occurrence of each nucleotide in cover_dna with encrypted_dna
#     result = list(cover_dna)
#     count = {}
#     enc_idx = 0
#     for i, n in enumerate(cover_dna):
#         count[n] = count.get(n, 0) + 1
#         if count[n] == 2 and enc_idx < len(encrypted_dna):
#             result[i] = encrypted_dna[enc_idx]
#             enc_idx += 1
#     return ''.join(result)

def find_adjacent_repeats_indices(dna_sequence):
    indices = []
    for i in range(1, len(dna_sequence)):
        if dna_sequence[i] == dna_sequence[i - 1]:
            indices.append(i)
    return indices


def hide_encrypted_message_in_dna(dna_cover, encrypted_msg, repeated_indices):
    """
    Replace the repeated characters in dna_cover at repeated_indices using encrypted_msg and Table 2 substitution rules.
    """
    dna_cover = list(dna_cover)  # convert to list for mutation
    for i, index in enumerate(repeated_indices):
        if i >= len(encrypted_msg):
            break  # no more message to hide
        msg_char = encrypted_msg[i]
        repeated_base = dna_cover[index]
        substitute = table2_substitution.get((msg_char, repeated_base), repeated_base)
        dna_cover[index] = substitute
    return ''.join(dna_cover)



def insert_delimited_into_dna(delimited_encrypted, dna_binary, key2):
    # print(f"Original DNA binary: {dna_binary}")
    # print(f"Delimited Encrypted Binary: {delimited_encrypted}")

    modified_binary = []
    encrypted_index = 0
    dna_index = 0

    while dna_index < len(dna_binary) and encrypted_index < len(delimited_encrypted):
        # Take two bits from DNA binary
        if dna_index < len(dna_binary) - 1:
            modified_binary.append(dna_binary[dna_index:dna_index + 2])
            dna_index += 2
        elif dna_index < len(dna_binary):
            modified_binary.append(dna_binary[dna_index])  # If last bit is left
            dna_index += 1

        # Insert one bit from the encrypted binary
        if encrypted_index < len(delimited_encrypted):
            modified_binary.append(delimited_encrypted[encrypted_index])
            encrypted_index += 1

    # Join the list into a final binary string
    modified_binary_str = ''.join(modified_binary)
    # print(f"Modified DNA binary after inserting encrypted result: {modified_binary_str}")

    return modified_binary_str

# Huffman Compression
class HuffmanNode:
    def __init__(self, char, freq):
        self.char = char
        self.freq = freq
        self.left = None
        self.right = None

    def __lt__(self, other):
        return self.freq < other.freq

def build_huffman_tree(text):
    frequency = defaultdict(int)
    for char in text:
        frequency[char] += 1
    
    priority_queue = [HuffmanNode(char, freq) for char, freq in frequency.items()]
    heapq.heapify(priority_queue)
    
    while len(priority_queue) > 1:
        left = heapq.heappop(priority_queue)
        right = heapq.heappop(priority_queue)
        merged = HuffmanNode(None, left.freq + right.freq)
        merged.left = left
        merged.right = right
        heapq.heappush(priority_queue, merged)
    
    return priority_queue[0]

def generate_huffman_codes(node, current_code="", codebook={}):
    if node is None:
        return
    
    if node.char is not None:
        codebook[node.char] = current_code
    
    generate_huffman_codes(node.left, current_code + "0", codebook)
    generate_huffman_codes(node.right, current_code + "1", codebook)
    
    return codebook

def huffman_compress(text):
    tree = build_huffman_tree(text)
    codebook = generate_huffman_codes(tree)
    compressed_data = ''.join(codebook[char] for char in text)
    return compressed_data, codebook

# Reed-Solomon Error Correction
def reed_solomon_encode(data, ecc_symbols=10):
    rsc = RSCodec(ecc_symbols)
    return rsc.encode(data)

# Hamming Error Correction
def hamming_encode(data):
    encoded_data = ""
    for byte in data:
        bits = format(byte, '08b')
        p1 = int(bits[0]) ^ int(bits[1]) ^ int(bits[3]) ^ int(bits[4]) ^ int(bits[6])
        p2 = int(bits[0]) ^ int(bits[2]) ^ int(bits[3]) ^ int(bits[5]) ^ int(bits[6])
        p4 = int(bits[1]) ^ int(bits[2]) ^ int(bits[3]) ^ int(bits[7])
        p8 = int(bits[4]) ^ int(bits[5]) ^ int(bits[6]) ^ int(bits[7])
        
        encoded_bits = f"{p1}{p2}{bits[0]}{p4}{bits[1]}{bits[2]}{bits[3]}{p8}{bits[4]}{bits[5]}{bits[6]}{bits[7]}"
        encoded_data += encoded_bits
    
    return encoded_data

# AES Encryption (CBC Mode)
def aes_encrypt_cbc(binary_text, key):
    binary_bytes = int(binary_text, 2).to_bytes((len(binary_text) + 7) // 8, byteorder='big')
    padded_data = pad(binary_bytes, AES.block_size)
    cipher = AES.new(key, AES.MODE_CBC)
    iv = cipher.iv
    cipher_text = cipher.encrypt(padded_data)
    return iv + cipher_text

# AES Encryption (GCM Mode)
def aes_encrypt_gcm(binary_text, key):
    binary_bytes = int(binary_text, 2).to_bytes((len(binary_text) + 7) // 8, byteorder='big')
    cipher = AES.new(key, AES.MODE_GCM)
    cipher_text, tag = cipher.encrypt_and_digest(binary_bytes)
    return cipher.nonce + tag + cipher_text

# RSA Key Generation and Signing
def generate_rsa_keys():
    private_key = rsa.generate_private_key(public_exponent=65537, key_size=2048)
    public_key = private_key.public_key()
    
    private_pem = private_key.private_bytes(
        encoding=serialization.Encoding.PEM,
        format=serialization.PrivateFormat.PKCS8,
        encryption_algorithm=serialization.NoEncryption()
    )
    
    public_pem = public_key.public_bytes(
        encoding=serialization.Encoding.OpenSSH,
        format=serialization.PublicFormat.OpenSSH
    )
    
    return private_key, public_key, private_pem.decode(), public_pem.decode()

def sign_with_rsa(private_key, hashed_data):
    return private_key.sign(
        hashed_data,
        padding.PSS(mgf=padding.MGF1(hashes.SHA256()), salt_length=padding.PSS.MAX_LENGTH),
        hashes.SHA256()
    )

# ECDSA Key Generation and Signing
def generate_ecdsa_keys():
    private_key = ec.generate_private_key(ec.SECP256R1())
    public_key = private_key.public_key()
    
    private_pem = private_key.private_bytes(
        encoding=serialization.Encoding.PEM,
        format=serialization.PrivateFormat.PKCS8,
        encryption_algorithm=serialization.NoEncryption()
    )
    
    public_pem = public_key.public_bytes(
        encoding=serialization.Encoding.PEM,
        format=serialization.PublicFormat.SubjectPublicKeyInfo
    )
    
    return private_key, public_key, private_pem.decode(), public_pem.decode()

def sign_with_ecdsa(private_key, hashed_data):
    return private_key.sign(hashed_data, ec.ECDSA(hashes.SHA256()))

def hash_encoded_data(encoded_dna):
    digest = hashes.Hash(hashes.SHA256())
    digest.update(encoded_dna.encode())
    return digest.finalize()

def save_file(path, data, binary=False):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    mode = 'wb' if binary else 'w'
    with open(path, mode) as file:
        if binary:
            file.write(data if isinstance(data, bytes) else data.encode())
        else:
            file.write(data if isinstance(data, str) else data.decode())

def hex_to_dna(hex_str):
    """Convert hexadecimal string to DNA sequence using 2-bit encoding"""
    dna = []
    for c in hex_str.lower():
        # Convert hex character to 4-bit binary
        binary = format(int(c, 16), '04b')
        # Split into two 2-bit pairs and map to nucleotides
        dna.append(DNA_ENCODING_TABLE[binary[:2]])
        dna.append(DNA_ENCODING_TABLE[binary[2:]])
    return ''.join(dna)

def main():

    with open('encode_config.json', 'r') as f:
        config = json.load(f)
    
    secret_message = read_file(config["message_file"])
    dna_sequence = read_file(config["dna_file"])
    method = config["method"]
    error_correction = config["error_correction"]
    cipher_mode = config["cipher_mode"]
    signature_scheme = config["signature"]
    output_path = config.get("output", "encoded_dna.txt")
    key = config.get("key", 2)
    dna_seq_2 = read_file(config["dna_file_2"]) if method == "method2" else None

    start_total = time.time()

    # Generate keys based on selected signature scheme
    start = time.time()
    if signature_scheme == "ecdsa":
        private_key, public_key, private_pem, public_pem = generate_ecdsa_keys()
    else:
        private_key, public_key, private_pem, public_pem = generate_rsa_keys()
    print(f"Key Generation Time: {time.time() - start:.4f} seconds")

    # Huffman Compression
    start = time.time()
    compressed_binary, codebook = huffman_compress(secret_message)
    compressed_length = len(compressed_binary)  # Store length of compressed binary
    print(f"Huffman Compression Time: {time.time() - start:.4f} seconds")

    # Save compressed_length to a secure file
    save_file("output/secure/compressed_length.txt", str(compressed_length))
    # print(compressed_binary)

    # AES Encryption
    aes_key = os.urandom(32)
    start = time.time()
    if cipher_mode == "cbc":
        encrypted_data = aes_encrypt_cbc(compressed_binary, aes_key)
    else:
        encrypted_data = aes_encrypt_gcm(compressed_binary, aes_key)
    print(f"AES Encryption Time: {time.time() - start:.4f} seconds")

    # print("encrypted_binary",encrypted_data)

    # Error Correction
    start = time.time()
    if error_correction == "reed-solomon":
        processed_data1 = reed_solomon_encode(encrypted_data)
        processed_data = ''.join(format(byte, '08b') for byte in processed_data1)
    else:
        processed_data = hamming_encode(encrypted_data)
    print(f"Error Correction Time: {time.time() - start:.4f} seconds")

    delimited_encrypted = add_delimiters(processed_data)


    start = time.time()
    if method=='method1':
        # DNA Encoding
        dna_binary = dna_to_binary(dna_sequence)
        modified_dna_binary = insert_delimited_into_dna(delimited_encrypted, dna_binary, key)
        if len(modified_dna_binary) % 2 != 0:
            modified_dna_binary += '0'  # Pad with '0' to make length even
        encoded_dna = binary_to_dna(modified_dna_binary)
    else:
        # DNA Codon Mapping
        dna_coded = binary_to_codons(delimited_encrypted)

        # DNA XOR Encryption
        encrypted_dna = xor_dna_sequences(dna_coded, dna_seq_2[:len(dna_coded)])

        repeated_nucleo_indexes=find_adjacent_repeats_indices(dna_sequence)
        # Embedding in Ś
        encoded_dna = hide_encrypted_message_in_dna(dna_sequence, encrypted_dna, repeated_nucleo_indexes)
    print(f"DNA Embedding Time: {time.time() - start:.4f} seconds")

    # Digital Signature
    start = time.time()
    hashed_data = hash_encoded_data(encoded_dna)
    if signature_scheme == "ecdsa":
        signature = sign_with_ecdsa(private_key, hashed_data)
    else:
        signature = sign_with_rsa(private_key, hashed_data)
    print(f"Digital Signature Time: {time.time() - start:.4f} seconds")

    #hashing the signature
    start = time.time()
    signature_hash= hashlib.sha256(signature).hexdigest()
    hex_hash = signature_hash[:10]  # Take first 10 hex characters
    dna_hash = hex_to_dna(hex_hash)
    if len(encoded_dna) > 10:
        encoded_dna = encoded_dna[:10] + dna_hash + encoded_dna[10:]
    print(f"Digital Signature hash Time: {time.time() - start:.4f} seconds")

    # Verify signature for demonstration
    # verify_signature(public_key, hashed_data, signature)

    total_time = time.time() - start_total
    print(f"\nProcess completed in {total_time:.4f} seconds")
    # print("codebook:",codebook)
    # print("aes_key",aes_key)

        # Setup output directories
    os.makedirs("output/shared", exist_ok=True)
    os.makedirs("output/secure", exist_ok=True)

    # Save outputs
    save_file("output/shared/encoded_dna.txt", encoded_dna)
    save_file("output/shared/signature.bin", signature.hex())
    save_file("output/shared/public_key.pem", public_pem)

    save_file("output/secure/private_key.pem", private_pem)
    save_file("output/secure/aes_key.bin", aes_key, binary=True)
    save_dic("output/secure/codebook.json",codebook)

    print("\n=== Process Summary ===")
    print(f"method used: {method.upper()}")
    print(f"Error Correction: {error_correction.upper()}")
    print(f"Encryption Mode: {cipher_mode.upper()}")
    print(f"Signature Scheme: {signature_scheme.upper()}")
    print(f"Encoded DNA saved to: output/shared/encoded_dna.txt")
    print(f"Signature saved to: output/shared/signature.bin")
    print(f"Public key saved to: output/shared/public_key.pem")
    print(f"Private key saved to: output/secure/private_key.pem")
    print(f"AES key saved to: output/secure/aes_key.bin")
    print(f"codebook is saved to: output/secure/codebook.json")
    # print("compressed binary:",compressed_binary)

if __name__ == "__main__":
    main()

    
