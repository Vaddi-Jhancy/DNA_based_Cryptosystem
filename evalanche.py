def read_file(file_path):
    with open(file_path, 'r') as file:
        data = file.read().strip()  # Read the file and remove extra spaces/newlines
    return data

def dna_to_binary(dna_sequence):
    """Convert DNA sequence to binary using a simple encoding."""
    DNA_DECODING_TABLE = {'A': '00', 'T': '01', 'G': '10', 'C': '11'}
    return ''.join(DNA_DECODING_TABLE[nucleotide] for nucleotide in dna_sequence)

def avalanche_effect(original_dna, modified_dna):
    """Calculate the avalanche effect between two DNA sequences."""
    original_binary = dna_to_binary(original_dna)
    modified_binary = dna_to_binary(modified_dna)
    
    # Ensure both sequences are of equal length
    min_length = min(len(original_binary), len(modified_binary))
    original_binary = original_binary[:min_length]
    modified_binary = modified_binary[:min_length]
    
    # Count bit differences
    bit_changes = sum(1 for i in range(min_length) if original_binary[i] != modified_binary[i])
    
    # Compute avalanche effect percentage
    avalanche_percentage = (bit_changes / min_length) * 100
    
    return avalanche_percentage

# Example Usage
original_dna = read_file("dna_.txt") # Example DNA sequence
modified_dna = read_file("second_encoded.txt")  # Modified sequence with slight change

avalanche_score = avalanche_effect(original_dna, modified_dna)
print(f"Avalanche Effect: {avalanche_score:.2f}%")
