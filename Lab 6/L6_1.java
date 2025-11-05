import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class L6_1 {
    static class DNAFragment implements Comparable<DNAFragment> {
        String sequence;
        int length;

        public DNAFragment(String sequence) {
            this.sequence = sequence;
            this.length = sequence.length();
        }

        // Sorting based on length (for gel migration simulation)
        @Override
        public int compareTo(DNAFragment other) {
            return Integer.compare(this.length, other.length);
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

    public static List<DNAFragment> generateRandomFragments(String dnaSequence, int numSamples) {
        List<DNAFragment> fragments = new ArrayList<>();
        Random random = new Random();

        for (int i = 0; i < numSamples; i++) {
            // Ensure fragment length between 100-3000 bases
            int fragmentLength = random.nextInt(2901) + 100; // 100 to 3000

            int startPosition = random.nextInt(Math.max(1, dnaSequence.length() - fragmentLength + 1));

            String fragment = dnaSequence.substring(startPosition, startPosition + fragmentLength);
            fragments.add(new DNAFragment(fragment));
        }

        return fragments;
    }

    public static void visualizeGelElectrophoresis(List<DNAFragment> fragments) {
        fragments.sort(null);

        int maxGelLength = 100;
        int maxFragmentLength = fragments.getLast().length;

        System.out.println("\n Gel Electrophoresis Simulation");
        System.out.println("Gel Migration Visualization (Shortest -> Longest):");

        for (DNAFragment fragment : fragments) {
            int migrationDistance = (int) ((double) fragment.length / maxFragmentLength * maxGelLength);

            // Create visualization bar
            String bar = String.format("%3d: ", fragment.length) +
                    "=".repeat(Math.max(0, migrationDistance));

            System.out.println(bar);
        }
    }

    public static void main(String[] args) {
        try {
            String fastaFilePath = "rna.fna";

            String dnaSequence = readFASTAFile(fastaFilePath);

            List<DNAFragment> fragments = generateRandomFragments(dnaSequence, 10);

            visualizeGelElectrophoresis(fragments);

        } catch (IOException e) {
            System.err.println("Error reading FASTA file: " + e.getMessage());
        }
    }
}

