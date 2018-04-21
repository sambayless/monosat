package monosat;

/**
 * Unchecked runtime exception thrown when attempting to access the assignment of
 * a literal, bitvector, or graph property that does not have an assignment.
 */
public class NoModelException extends RuntimeException{
    public NoModelException(String message) {
        super(message);
    }
}