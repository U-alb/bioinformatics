public class Lab3_1 {

    public static void main(String[] args) {
        String DNA = "AATTCCGGGATC";
        int freqA = countBase(DNA, 'A');
        int freqT = countBase(DNA, 'T');
        int freqC = countBase(DNA, 'C');
        int freqG = countBase(DNA, 'G');

        double tms = calculateMeltingTemperatureSimple(freqA, freqT, freqC, freqG);
        double tmc = calculateMeltingTemperatureComplex(DNA, freqA, freqT, freqC, freqG);

        System.out.printf("Simple melting Temperature (Tm): %.2f°C%n", tms);
        System.out.printf("Complex Melting Temperature (Tm): %.2f°C%n", tmc);
    }

    static int countBase(String dna, char base) {
        int count = 0;
        for (char c : dna.toCharArray()) {
            if (Character.toUpperCase(c) == Character.toUpperCase(base)) {
                count++;
            }
        }
        return count;
    }

    static double calculateMeltingTemperatureSimple(int freqA, int freqT, int freqC, int freqG) {
        return 4 * (freqG + freqC) + 2 * (freqA + freqT);
    }

    static double calculateMeltingTemperatureComplex(String dna, int freqA, int freqT, int freqC, int freqG) {
        double f = (Math.log10(0.05));
        return 81.5 + 16.6 * f + 0.41 * ((double) (freqG + freqC) / dna.length())*100 - (double) (600 / dna.length());
    }
}
