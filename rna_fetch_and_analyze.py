# Imports
# Have Byopython and Python 3 installed!
from Bio import Entrez, SeqIO
from Bio.Seq import Seq

def fetch_sequence(accession_id):
    """
    Fetch nucleotide sequence from NCBI using the accession ID.
    """
    # Make sure to put your email. Use institutional email if possible.
    Entrez.email = "your_email@example.com"  # Replace with your email
    handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    return record

def calculate_gc_content(seq):
    """
    Calculate the GC content of an RNA sequence.
    """
    gc_count = seq.count('G') + seq.count('C')
    gc_content = (gc_count / len(seq)) * 100
    return gc_content

def get_reverse_complement(seq):
    """
    Get the reverse complement of an RNA sequence.
    """
    return seq.reverse_complement()

def convert_dna_to_rna(dna_seq):
    """
    Convert a DNA sequence to an RNA sequence by replacing 'T' with 'U'.
    """
    return dna_seq.replace('T', 'U').replace('t', 'u')

def save_to_file(filename, content):
    """
    Save the given content to a text file.
    """
    with open(filename, 'w') as file:
        file.write(content)

def main():
    # Input accession ID
    accession_id = input("Enter the accession ID of the nucleotide sequence: ")
    
    # Fetch the sequence from NCBI
    record = fetch_sequence(accession_id)
    sequence = str(record.seq)
    
    # Detect if the sequence is DNA and convert to RNA if necessary
    conversion_note = ""
    if 'T' in sequence and 'U' not in sequence:
        conversion_note = "The sequence was detected as DNA and converted to RNA.\n"
        sequence = convert_dna_to_rna(sequence)
    
    # Create a Seq object
    seq = Seq(sequence)
    
    # Calculate GC content
    gc_content = calculate_gc_content(seq)
    
    # Get reverse complement and convert to RNA
    reverse_complement = get_reverse_complement(seq).transcribe()
    
    # Prepare content for output file
    output_content = (
        f"Fetched sequence: {record.description}\n\n"
        f"Original Sequence:\n{sequence}\n\n"
        f"GC Content: {gc_content:.2f}%\n\n"
        f"Reverse Complement:\n{reverse_complement}\n\n"
        f"{conversion_note}"
    )
    
    # Create dynamic output filename
    output_filename = f"sequence_analysis_{accession_id}.txt"
    save_to_file(output_filename, output_content)
    
    # Print the description and GC content
    print(f"Fetched sequence: {record.description}")
    print(f"GC Content: {gc_content:.2f}%")

if __name__ == "__main__":
    main()
