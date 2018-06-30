package monosat;

import java.util.ArrayList;


/**
 * Produces MonoSAT constraints for testing purposes.
 * These are not intended to be demonstrations of how to write constraints
 * effectively in the solver.
 */
public class Constraints {
    public static ArrayList<ArrayList<Lit>> nqueens(Solver s, int n){
        ArrayList<ArrayList<Lit>> rows = new ArrayList<ArrayList<Lit>>();
        for(int i = 0;i<n;i++){
            ArrayList<Lit> row = new ArrayList<Lit>();
            for(int j = 0;j<n;j++){
                row.add(new Lit(s));
            }
            rows.add(row);
        }

        ArrayList<ArrayList<Lit>> cols = new ArrayList<ArrayList<Lit>>();
        for(int i = 0;i<n;i++){
            ArrayList<Lit> col = new ArrayList<Lit>();
            for(int j = 0;j<n;j++){
                col.add(rows.get(j).get(i));
            }
            cols.add(col);
        }



        for(int i = 0;i<n;i++){
            ArrayList<Lit> row = rows.get(i);
            //each row must have a queen
            s.addClause(row);
            //each row must have exactly one queen
            s.assertAtMostOne(row);
        }

        for(int i = 0;i<n;i++){
            ArrayList<Lit> col = cols.get(i);
            //each col must have a queen
            s.addClause(col);
            //each col must have exactly one queen
            s.assertAtMostOne(col);
        }

        //each diagonal must have at most one queen
        for(int i = 0;i<n;i++){
            for(int j = 0;j<n;j++){
                for(int k = 1;k<n-i && k<n-j;k++){
                    assert(i+k<n);
                    assert(j+k<n);
                    s.assertAtMostOne(rows.get(i).get(j), rows.get(i+k).get(j+k));
                }
                for(int k = 1;k<=i && k<=j;k++){
                    s.assertAtMostOne(rows.get(i).get(j), rows.get(i-k).get(j-k));
                }
                for(int k = 1;k<= i && k< n-j;k++){
                    s.assertAtMostOne(rows.get(i).get(j), rows.get(i-k).get(j+k));
                }
                for(int k = 1;k<= j && k< n-i;k++){
                    s.assertAtMostOne(rows.get(i).get(j), rows.get(i+k).get(j-k));
                }
            }
        }
        return rows;
    }

    public static void unsatQueens(Solver s, int n){
        ArrayList<ArrayList<Lit>> rows = nqueens(s,n);
        ArrayList<Lit> lits = new ArrayList<>();
        rows.forEach(lits::addAll);

        s.assertPB(lits,Comparison.GEQ,n+1);//assert n+1 queens are enabled
    }
}
