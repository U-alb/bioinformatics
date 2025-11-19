from Bio import SeqIO

def find_inverted_repeats(sequence, min_length=4, max_length=6):
    repeats_info = []
    for length in range(min_length, max_length + 1):
        for i in range(len(sequence) - length + 1):
            subseq = sequence[i:i + length]
            inverted_subseq = subseq[::-1]
            # Check for inverted repeat in the sequence
            pos = sequence.find(inverted_subseq, i + length)
            if pos != -1:
                repeats_info.append((subseq, i, pos + length - 1))  # Save (inverted_repeat, start, end)
    return repeats_info

fasta_file = "sequence.fasta"
inverted_repeats_report = {}

for record in SeqIO.parse(fasta_file, "fasta"):
    genome_id = record.id
    sequence = str(record.seq)
    inverted_repeats_report[genome_id] = find_inverted_repeats(sequence)

with open("inverted_repeats_report.txt", "w") as report:
    for genome_id, repeats in inverted_repeats_report.items():
        report.write(f"Genome ID: {genome_id}\n")
        number_of_repeats = len(repeats)
        report.write(f"Total Transposons Found: {number_of_repeats}\n")

        if repeats:
            for repeat, start, end in repeats:
                report.write(f"Inverted Repeat: {repeat} | Start: {start} | End: {end}\n")
        else:
            report.write("No inverted repeats found\n")

        report.write("\n")