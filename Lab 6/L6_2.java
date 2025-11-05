import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class L6_2 {
    static class RestrictionEnzyme {
        String name;
        String recognitionSite;

        public RestrictionEnzyme(String name, String recognitionSite) {
            this.name = name;
            this.recognitionSite = recognitionSite;
        }

        public List<String> digest(String dnaSequence) {
            List<String> fragments = new ArrayList<>();
            Pattern pattern = Pattern.compile(recognitionSite);
            Matcher matcher = pattern.matcher(dnaSequence);

            int lastEnd = 0;
            while (matcher.find()) {
                // Extract fragment between cuts
                String fragment = dnaSequence.substring(lastEnd, matcher.start());
                if (!fragment.isEmpty()) {
                    fragments.add(fragment);
                }
                lastEnd = matcher.end();
            }

            // Add the final fragment
            String finalFragment = dnaSequence.substring(lastEnd);
            if (!finalFragment.isEmpty()) {
                fragments.add(finalFragment);
            }

            return fragments;
        }
    }

    public static String readFASTAFile(String filePath) throws IOException {
        StringBuilder sequence = new StringBuilder();
        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line;
            boolean isSequenceLine = false;
            while ((line = br.readLine()) != null) {
                if (line.startsWith(">")) {
                    isSequenceLine = true;
                    continue;
                }
                if (isSequenceLine) {
                    sequence.append(line.trim());
                }
            }
        }
        return sequence.toString().toUpperCase();
    }

    public static void visualizeRestrictionDigestion(List<RestrictionEnzyme> enzymes, String dnaSequence) {
        System.out.println("\nRestriction Enzyme Digestion Simulation");

        // Table header
        System.out.printf("%-15s | %-20s | %-10s | %s\n",
                "Enzyme", "Recognition Site", "Fragments", "Fragment Lengths");
        System.out.println("-".repeat(70));

        for (RestrictionEnzyme enzyme : enzymes) {
            // Digest the sequence with current enzyme
            List<String> fragments = enzyme.digest(dnaSequence);

            // Prepare fragment lengths
            List<Integer> fragmentLengths = fragments.stream()
                    .map(String::length)
                    .toList();

            // Print results
            System.out.printf("%-15s | %-20s | %-10d | %s\n",
                    enzyme.name,
                    enzyme.recognitionSite,
                    fragments.size(),
                    fragmentLengths
            );
        }
    }

    public static void main(String[] args) {
        try {
            String fastaFilePath = "rna.fna";
            String dnaSequence = readFASTAFile(fastaFilePath);

            // Define restriction enzymes
            List<RestrictionEnzyme> enzymes = List.of(
                    new RestrictionEnzyme("EcoRI", "G[AATtc]TC"),      // Recognizes G^AATTC
                    new RestrictionEnzyme("BamHI", "G[GATCC]C"),        // Recognizes G^GATCC
                    new RestrictionEnzyme("HindIII", "A[AGCTT]T"),     // Recognizes A^AGCTT
                    new RestrictionEnzyme("XhoI", "C[TCGAG]G"),        // Recognizes C^TCGAG
                    new RestrictionEnzyme("NotI", "GC[GGCCGC]C")       // Recognizes GC^GGCCGC
            );

            // Visualize restriction digestion
            visualizeRestrictionDigestion(enzymes, dnaSequence);

        } catch (IOException e) {
            System.err.println("Error reading FASTA file: " + e.getMessage());
        }
    }
}
