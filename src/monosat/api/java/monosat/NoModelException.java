package monosat;

import java.util.NoSuchElementException;
/**
 * Unchecked runtime exception thrown when attempting to access the assignment of
 * a literal, bitvector, or graph property that does not have an assignment.
 */
public class NoModelException extends NoSuchElementException {
    public NoModelException(String message) {
        super(message);
    }
}