# DNA_based_Cryptosystem

DNA‑Based Secure Data Hiding System
This project implements a DNA‑based secure data hiding system that combines compression, encryption, error correction, steganography, and digital signatures to provide confidentiality, integrity, authenticity, and high robustness for secret messaging using DNA sequences as the cover medium.​

**Project overview**

The system hides a compressed and encrypted secret message inside a DNA sequence so that the modified DNA still appears like a normal biological sequence, making detection extremely difficult while also protecting the content cryptographically. Two DNA‑steganography methods are provided, along with strong AES‑256 encryption and RSA/ECDSA signatures to give a multilayered security architecture suitable for secure messaging and DNA data storage.​

**Key features**

*DNA steganography: Secret bits are embedded into organism DNA so that even large payloads are statistically and visually indistinguishable from normal DNA, especially when long sequences (up to billions of bases) are used.​

*Strong encryption (AES‑256, CBC/GCM): The message is encrypted with AES‑256 in CBC or GCM mode, giving a key‑space of 2^256 possibilities and ensuring unreadability without the symmetric key, even if the embedded data is extracted.​

*Error correction (Reed–Solomon, Hamming): Additional parity/check symbols are added to detect and correct both single‑bit and burst errors that may arise during DNA storage, transmission, or sequencing.​

*Digital signatures (RSA/ECDSA): SHA‑256 hashing plus RSA‑2048 or ECDSA signatures provide sender authentication and integrity verification before decryption is attempted.​

*High capacity and low modification rate: Human DNA scale (up to 1.6 billion bases) allows very high hiding capacity with extremely low relative modification of the original sequence, preserving biological functionality.​

**Methodologies**

1. Compression (Huffman coding)
The plaintext message is first compressed using lossless Huffman coding to reduce size by roughly 30–50%, mapping frequent characters to shorter bit codes and storing the codebook separately for decompression. This step both increases capacity and reduces the number of bits that must be encrypted, corrected, and embedded.​

2. Encryption (AES‑256 in CBC/GCM)
AES‑256 is applied on the compressed bitstream in 256‑bit blocks using a secret symmetric key.​

CBC mode: Uses a 256‑bit IV; each plaintext block is XORed with the previous ciphertext (or IV for the first block) before encryption, producing only ciphertext output.​

GCM mode: Uses a 96‑bit nonce and a counter to generate a keystream, encrypting in a stream‑like fashion and producing both ciphertext and an authentication tag for integrity.​

3. Error correction (Reed–Solomon, Hamming)
Reed–Solomon codes add check symbols (typically at the beginning) that can correct burst errors over blocks of bits.​

Hamming codes add 4 parity bits per 8‑bit byte (12‑bit codewords) to correct single‑bit errors at the cost of increased length.​
These codes are removed after decoding and correction during decryption.​

**DNA encoding and embedding**

Before embedding, fixed delimiters of 7 bytes are added to the start and end of the protected bitstream so that the payload can be located inside the DNA during extraction. DNA encoding uses a simple 2‑bit mapping such as A = 00, T = 01, G = 10, C = 11 to convert between bits and bases.​

**Method 1 – Insertion‑based embedding**

Bits (with delimiters and ECC) are inserted into the DNA bitstream at positions determined by an embedding key (e.g., every k bits), producing a longer binary string.​
The modified bitstring is then mapped back to bases, yielding an encoded DNA sequence whose length is larger than the original (e.g., more than 15,000 bases in experiments).​

Strengths:​

Single‑layer security (one cover DNA sequence).

Fully blind: receiver does not need the original DNA, only the encoded DNA and keys.

Higher payload capacity and usually simpler/faster to implement.

**Method 2 – Substitution and XOR‑based embedding**

The protected bits are grouped into 6‑bit chunks and mapped to codons (3 bases) using a codon mapping table.​
Codons are decoded back to bits and XORed with a second DNA sequence’s bits; the result is re‑encoded into bases and then substituted into the main (human) DNA using a substitution table, typically on repeated bases.​

The overall DNA length remains the same (e.g., 15,000 bases), as bases are replaced, not inserted.​

Strengths:​

Double‑layer security (cover DNA + second DNA).

Partially blind: receiver needs additional information (indices, substitution table, second DNA).

Lower capacity and higher complexity, but stronger confidentiality and confounding structure.

**Digital signatures and authenticity**

A SHA‑256 hash of the final encoded DNA is generated.​

This hash is signed using either RSA‑2048 or ECDSA to create a digital signature, where RSA uses large primes and yields a ~2048‑bit signature, while ECDSA achieves similar security with ~512‑bit signatures and better performance.​

The signature (or its hash) is itself converted to DNA bases and embedded, enabling verification that both the sender and the DNA sequence are genuine and untampered.​

*End‑to‑end workflow*

Encryption / embedding
Compress plaintext with Huffman coding and store the character–code mapping (codebook).​

Encrypt the compressed data using AES‑256 in CBC or GCM mode, obtaining ciphertext (and tag in GCM).​

Add parity/check bits using Reed–Solomon (for burst errors) or Hamming (for single‑bit errors).​

Add start and end delimiters of 7 bytes to frame the protected payload.​

Embed the resulting bits into DNA using Method 1 or Method 2 to obtain the encoded DNA.​

Compute SHA‑256 of the encoded DNA, sign with RSA/ECDSA private key, convert signature to DNA, and embed it as well.​

Public information: encoded DNA and public key(s); private information: AES key, codebook, embedding key(s), and signing private key.​

Decryption / extraction
Extract and decode the embedded signature, then verify it using the sender’s public key; abort if verification fails.​

Reverse Method 1 or Method 2 to recover the raw embedded bitstream from encoded DNA.​

Locate delimiters and extract only the protected payload between them.​

Apply the appropriate error‑correction decoding (Reed–Solomon or Hamming) to fix any bit errors.​

Decrypt the corrected ciphertext with AES‑256 in the correct mode (CBC/GCM) using the shared key.​

Decompress the decrypted data using the Huffman codebook to reconstruct the original message.​

*Experimental setup and results*
Example cover DNA length: 15,000 bases (30,000 bits with 2‑bit encoding), with the conceptual design targeting sequences up to 1.6 billion bases.​

Example secret message length: 614 ASCII characters (4,912 bits), compressed to ≈2,560 bits (~50% reduction) using Huffman coding.​

After encryption, error correction, IV/nonce, tag, padding, and delimiters, the final protected payload is ≈3,000 bits before embedding.​

Method 1 increases DNA length beyond 15,000 bases, while Method 2 preserves the original DNA length.​

The system demonstrates correct recovery of the original message under both single‑bit and burst error scenarios (Hamming vs Reed–Solomon) and maintains very low modification rates when large DNA sequences are used.​

Security and evaluation
Key evaluation factors include cracking probability, layers of security, blindness, modification rate, payload capacity, preservation of DNA functionality, hiding capacity, speed, confidentiality, authentication, integrity, correctness, and avalanche effect.​

Cracking probability: AES‑256 encryption gives a key‑space of 2^256, making brute‑force infeasible; Method 2’s multi‑dependency structure further reduces practical cracking probability.​

Security layers: Method 1 uses a single cover DNA, while Method 2 uses a cover plus second DNA, providing double‑layer obfuscation.​

Blindness: Method 1 is fully blind (no original DNA required at the receiver), whereas Method 2 is partially blind due to dependencies on auxiliary data.​

Modification rate and DNA functionality: Only a small fraction of bases are modified; with very long DNA sequences, the functional impact remains negligible, preserving protein‑coding roles.​

Capacity: Method 1 offers higher payload and nearly unlimited capacity in very long sequences; Method 2 is more constrained by the distribution of repeated bases.​

Correctness and robustness: Hamming codes ensure single‑bit correction, while Reed–Solomon handles burst errors, allowing reliable message recovery in realistic noisy conditions.​

Avalanche effect: The ratio of DNA bits to message bits keeps observable avalanche effects below about 5%, making small input changes still lead to sufficiently unpredictable outputs in the steganographic context.​

Use cases
This DNA‑based secure data hiding system is suitable for:​

Long‑term archival of sensitive data within synthetic or biological DNA storage systems.​

High‑confidentiality messaging where both secrecy and undetectability are required.​

Scenarios demanding authenticity and tamper‑evidence through integrated digital signatures.​

The design targets a balance between practicality (through Method 1), stronger multi‑layer security (through Method 2), and robustness via standard cryptographic and coding primitives, making it a comprehensive framework for secure DNA steganography.
