package io.github.twalgor.main;

import io.github.twalgor.common.*;
import io.github.twalgor.greedy.MMAF;
import io.github.twalgor.recursive.ContractionRecursive;
import io.github.twalgor.safesep.RootedMinorBoundedDepthBacktrack;
import io.github.twalgor.safesep.SafeSepDecomposition;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Date;

public class Main {
  static final int VERSION = 2;
  
  private static void solve(String name, String graphPath, String certPath) {
    Calendar cl = Calendar.getInstance();
    Date date = new Date();
    
    date = cl.getTime();
    System.out.println(date);
    
    File graphDir = new File(graphPath);
    File graphFile = new File(graphDir, name + ".gr");

    Graph g = Graph.readGraph(graphFile);
    System.out.println(name + " n " + g.n);


    System.out.println("Testing if input is okay");
    g.printAdjacencyList();

    long t0 = System.currentTimeMillis();
    SafeSepDecomposition ssd = new SafeSepDecomposition(g);
    System.out.println("Initialized SafeSepDecomposition object");
    TreeDecomposition std = ssd.decomposeSafe();

    std.printDecomposition();
    
    Integer[] bx = new Integer[std.nb];
    for (int b = 1; b <= std.nb; b++) {
      bx[b - 1] = b;
    }
    
    Arrays.sort(bx, (b1, b2) -> -(std.bags[b1].length - std.bags[b2].length));

    System.out.println("Sorted the bags of the tree decomposition by size (descending)");
    Graph tr = g.copy();
    int twMax = 0;
    int bxMax = 0;
    LocalGraph lgMax = null;
    ContractionRecursive crMax = null;
    System.out.println("Now iterating through the sorted bags");


    for (int i = 0; i < bx.length; i++) {

      System.out.println("Processing bag " + bx[i] + " with size " + std.bags[bx[i]].length);
      LocalGraph lg = new LocalGraph(g, std.bagAt(bx[i]));
      lg.printVertices();
      if (i > 0)
      {
        Graph f = lg.h.copy();
        MMAF mmaf = new MMAF(f);
        mmaf.triangulate();
        System.out.println("After triangulation the width is: " + mmaf.width);
        if (mmaf.width <= twMax)
        {
          System.out.println("MMAF triangulation width is less than maximum Treewidth: "+ twMax);
          Chordal chordal = new Chordal(f);

          int cliqueIndex = 1;

          for (XBitSet clique: chordal.maximalCliques()) {
            System.out.println("Clique Number: "+ cliqueIndex);

            System.out.print("Clique " + cliqueIndex + ": ");
            for (int vertex = clique.nextSetBit(0); vertex >= 0; vertex = clique.nextSetBit(vertex + 1)) {
              System.out.print(vertex + " ");
            }
            System.out.println();
            cliqueIndex++;
            tr.fill(clique.convert(lg.inv));
          }
          continue;
        }
      }
      ContractionRecursive cr = new ContractionRecursive(lg.h);
      System.out.println("Contraction Recursive Called on the subgraph");
      TreeDecomposition td = cr.decompose();
      td.printDecomposition();

      for (XBitSet bag: td.setOfBags()) {
        XBitSet bag1 = bag.convert(lg.inv);
        assert g.isPMC(bag1);
        tr.fill(bag1);
      }
      if (td.width > twMax) {
        twMax = td.width;
        bxMax = bx[i];
        crMax = cr;
        lgMax = lg;
      }
    }

    System.out.println("Converted the graph to a tree decomposition using Chordal.chordalToTD");
    TreeDecomposition td = Chordal.chordalToTD(tr);
    td.g = g;
    
    Minor minor0 = localizingMinor(std, bxMax, lgMax);
    System.out.println("Computed and combined the localizing minor with recursive contraction results");
    Minor minor = crMax.obs.composeWith(minor0);
    long t = System.currentTimeMillis();
    System.out.println("certificates computed " + (t - t0) + " millisecs");
    
    td.validate();
    
    File certDir =  new File(certPath);
    certDir.mkdirs();
    
    File certFile =  new File(certDir, name + ".twc");
    
    try {
      PrintStream ps = new PrintStream(new FileOutputStream(certFile, false));
      ps.println("c treewidth certificate " + name + 
          " n " + g.n + " m " + g.numberOfEdges() + " graphpath " + graphPath);
      printParams(ps);
      ps.println("c width " + td.width + " time(ms) " + (t - t0));
      ps.println("s nbags " + td.nb);
      for (int b = 1; b <= td.nb; b++) {
        ps.print("b " + b);
        for (int i: td.bags[b]) {
          ps.print(" " + (i + 1));
        }
        ps.println();
      }
      for (Edge e: td.toTree().edgeList()) {
        ps.println(e.u + " " + e.v);
      }
      ps.println("s ncomponents " + minor.components.length);
      for (int i = 0; i < minor.components.length; i++) {
        ps.print((i + 1));
        for (int v: minor.components[i].toArray()) {
          ps.print(" " + (v + 1));
        }
        ps.println();
      }
      
    } catch (FileNotFoundException e) {
      e.printStackTrace();
    }

  }

  static void printParams(PrintStream ps) {
    ps.print("c ");
    ps.print(" UNIT_BUDGET_UB " + ContractionRecursive.UNIT_BUDGET_UB);
    ps.print(" N_CONDUCIVES_RATIO " + ContractionRecursive.N_CONDUCIVES_RATIO);
    ps.println();
  }

  static Minor localizingMinor(TreeDecomposition std, int b, LocalGraph lg) {
    XBitSet[] ca = new XBitSet[lg.h.n];
    for (int i = 0; i < ca.length; i++) {
      ca[i] = new XBitSet(new int[] {lg.inv[i]});
    }
    
    Graph g = std.g;
    XBitSet bag = std.bagAt(b);
    ArrayList<XBitSet> components = g.separatedComponents(bag);
    for (XBitSet compo: components) {
      XBitSet sep = g.neighborSet(compo);
      Subgraph sub = new Subgraph(g, sep.unionWith(compo));
      XBitSet root = sep.convert(sub.conv);
      RootedMinorBoundedDepthBacktrack rmbdb = 
          new RootedMinorBoundedDepthBacktrack(sub.h, root);
      boolean cliqueMinorFound = rmbdb.hasCliqueMinor();
      assert cliqueMinorFound;
      Minor rMinor = rmbdb.rMinor;
      for (XBitSet c: rMinor.components) {
        XBitSet intersect = c.intersectWith(root);
        assert intersect.cardinality() == 1;
        int v = intersect.nextSetBit(0);
        int w = sub.inv[v];
        ca[lg.conv[w]].or(c.convert(sub.inv));
      }
    }
        
    return new Minor(g, ca);
  }

  static String vaToString(int[] va) {
    StringBuilder sb = new StringBuilder();
    for (int v: va) {
      sb.append(" " + (v + 1));
    }
    return sb.toString();
  }
  
  static String baToString(int b, int[] ba) {
    StringBuilder sb = new StringBuilder();
    return sb.toString();
  }
    
  public static void main(String[] args) {
    if (args.length == 0) {
//      String name = "Promedas_68_13";
//      String name = "Promedas_28_10";
//      String name = "Promedas_11_7";
//      String name = "Promedas_27_8";
      //String name = "Promedus_27_15";
//      String name = "Pedigree_11_6";
//      String name = "Pedigree_12_12";
//      String name = "Promedus_38_15";
//      String name = "Promedus_38_14";
//      String name = "MD5-32-1.gaifman_4";
      String name = "ex001";
      args = new String[] {
          name, 
          "instance\\PACE2017bonus_gr",
          "twcert\\PACE2017bonus_gr"
      };
    }
    assert args.length >= 3;
    solve(args[0], args[1], args[2]);
  }
}
