import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

public class Lab4_1 {
    private static final Map<String, String> geneticCode = new HashMap<>();

    static {
        geneticCode.put("UUU", "Phe");
        geneticCode.put("UUC", "Phe");
        geneticCode.put("UUA", "Leu");
        geneticCode.put("UUG", "Leu");
        geneticCode.put("UCU", "Ser");
        geneticCode.put("UCC", "Ser");
        geneticCode.put("UCA", "Ser");
        geneticCode.put("UCG", "Ser");
        geneticCode.put("UAU", "Tyr");
        geneticCode.put("UAC", "Tyr");
        geneticCode.put("UAA", "Stop");
        geneticCode.put("UAG", "Stop");
        geneticCode.put("UGU", "Cys");
        geneticCode.put("UGC", "Cys");
        geneticCode.put("UGA", "Stop");
        geneticCode.put("UGG", "Trp");
        geneticCode.put("CUU", "Leu");
        geneticCode.put("CUC", "Leu");
        geneticCode.put("CUA", "Leu");
        geneticCode.put("CUG", "Leu");
        geneticCode.put("CCU", "Pro");
        geneticCode.put("CCC", "Pro");
        geneticCode.put("CCA", "Pro");
        geneticCode.put("CCG", "Pro");
        geneticCode.put("CAU", "His");
        geneticCode.put("CAC", "His");
        geneticCode.put("CAA", "Gln");
        geneticCode.put("CAG", "Gln");
        geneticCode.put("CGU", "Arg");
        geneticCode.put("CGC", "Arg");
        geneticCode.put("CGA", "Arg");
        geneticCode.put("CGG", "Arg");
        geneticCode.put("AUU", "Ile");
        geneticCode.put("AUC", "Ile");
        geneticCode.put("AUA", "Ile");
        geneticCode.put("AUG", "Met");
        geneticCode.put("ACU", "Thr");
        geneticCode.put("ACC", "Thr");
        geneticCode.put("ACA", "Thr");
        geneticCode.put("ACG", "Thr");
        geneticCode.put("AAU", "Asn");
        geneticCode.put("AAC", "Asn");
        geneticCode.put("AAA", "Lys");
        geneticCode.put("AAG", "Lys");
        geneticCode.put("AGU", "Ser");
        geneticCode.put("AGC", "Ser");
        geneticCode.put("AGA", "Arg");
        geneticCode.put("AGG", "Arg");
        geneticCode.put("GUU", "Val");
        geneticCode.put("GUC", "Val");
        geneticCode.put("GUA", "Val");
        geneticCode.put("GUG", "Val");
        geneticCode.put("GCU", "Ala");
        geneticCode.put("GCC", "Ala");
        geneticCode.put("GCA", "Ala");
        geneticCode.put("GCG", "Ala");
        geneticCode.put("GAU", "Asp");
        geneticCode.put("GAC", "Asp");
        geneticCode.put("GAA", "Glu");
        geneticCode.put("GAG", "Glu");
        geneticCode.put("GGU", "Gly");
        geneticCode.put("GGC", "Gly");
        geneticCode.put("GGA", "Gly");
        geneticCode.put("GGG", "Gly");
    }

    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);
        System.out.println("Enter the coding region of the gene (DNA sequence):");
        String dnaSequence = scanner.nextLine().toUpperCase();

        String rnaSequence = transcribeToRNA(dnaSequence);
        String aminoAcidSequence = translateToAminoAcids(rnaSequence);

        System.out.println("RNA Sequence: " + rnaSequence);
        System.out.println("Amino Acid Sequence: " + aminoAcidSequence);
    }

    private static String transcribeToRNA(String dna) {
        return dna.replace('T', 'U');
    }

    private static String translateToAminoAcids(String rna) {
        StringBuilder aminoAcidSequence = new StringBuilder();

        for (int i = 0; i < rna.length(); i += 3) {
            if (i + 2 < rna.length()) {
                String codon = rna.substring(i, i + 3);
                String aminoAcid = geneticCode.get(codon);
                if (aminoAcid != null) {
                    if (aminoAcid.equals("Stop")) {
                        break;
                    }
                    aminoAcidSequence.append(aminoAcid).append("-");
                }
            }
        }

        if (!aminoAcidSequence.isEmpty()) {
            aminoAcidSequence.setLength(aminoAcidSequence.length() - 1);
        }

        return aminoAcidSequence.toString();
    }
}