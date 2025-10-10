import java.util.*;

public class Lab2_1 {
    private static final char[] BASES = {'A', 'C', 'G', 'T'};

    static void main(String[] args) {
        String S = "ATTGTCCCAATCTGTTG".toUpperCase();

        Map<String, Integer> diCounts = initAll(2);
        Map<String, Integer> triCounts = initAll(3);

        int len = S.length();
        int diWindows = Math.max(0, len - 1);
        for (int i = 0; i + 2 <= len; i++) {
            String sub = S.substring(i, i + 2);
            if (isValid(sub)) diCounts.put(sub, diCounts.get(sub) + 1);
        }
        int triWindows = Math.max(0, len - 2);
        for (int i = 0; i + 3 <= len; i++) {
            String sub = S.substring(i, i + 3);
            if (isValid(sub)) triCounts.put(sub, triCounts.get(sub) + 1);
        }

        System.out.println("Sequence: " + S);
        System.out.println();

        System.out.println("Dinucleotide count and frequencies (count / " + diWindows + "):");
        for (String k : sortedKeys(diCounts.keySet())) {
            int c = diCounts.get(k);
            double freq = diWindows > 0 ? (double) c / diWindows : 0.0;
            System.out.printf("%s\t%d\t%.6f%n", k, c, freq);
        }

        System.out.println();
        System.out.println("Trinucleotide count and frequencies (count / " + triWindows + "):");
        for (String k : sortedKeys(triCounts.keySet())) {
            int c = triCounts.get(k);
            double freq = triWindows > 0 ? (double) c / triWindows : 0.0;
            System.out.printf("%s\t%d\t%.6f%n", k, c, freq);
        }
    }

    private static Map<String, Integer> initAll(int width) {
        Map<String, Integer> map = new LinkedHashMap<>();
        generateRecursive(new StringBuilder(), width, map);
        return map;
    }

    private static void generateRecursive(StringBuilder sb, int width, Map<String, Integer> map) {
        if (sb.length() == width) {
            map.put(sb.toString(), 0);
            return;
        }
        for (char b : BASES) {
            sb.append(b);
            generateRecursive(sb, width, map);
            sb.setLength(sb.length() - 1);
        }
    }

    private static boolean isValid(String s) {
        for (int i = 0; i < s.length(); i++) {
            char c = s.charAt(i);
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T') return false;
        }
        return true;
    }

    private static List<String> sortedKeys(Set<String> keys) {
        List<String> list = new ArrayList<>(keys);
        Collections.sort(list);
        return list;
    }
}
