package io.github.twalgor.common;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public class LocalGraph {
  protected Graph g;
  public XBitSet vertices;
  public int[] conv;
  public int[] inv;
  public Graph h;
  public ArrayList<XBitSet> cliques;
  public Map<XBitSet, XBitSet> cliqueMap;
  
  public LocalGraph(Graph g, XBitSet vertices) {
    this(g, vertices, null);
  }
  
  public LocalGraph(Graph g, XBitSet vertices, ArrayList<XBitSet> givenCliques) {
    this.g = g;
    this.vertices = vertices;
    
    conv = new int[g.n];
    inv = new int[vertices.cardinality()];
    Arrays.fill(conv, -1);
    int w = 0;
    for (int v = vertices.nextSetBit(0); v >= 0; v = vertices.nextSetBit(v + 1)) {
      inv[w] = v;
      conv[v] = w++;
    }

    h = new Graph(w);
    h.inheritEdges(g, conv, inv);
    

    cliques = new ArrayList<> ();
    cliqueMap = new HashMap<>();

    if (givenCliques != null) {
      for (XBitSet givenClique: givenCliques) {
        XBitSet clique = givenClique.convert(conv);
        cliques.add(clique);
        h.fill(clique);
      }
    }
    else {
      ArrayList<XBitSet> components = g.separatedComponents(vertices);
//      System.out.println(components.size() + " components");
//      System.out.println("vertices: " + vertices);
//      System.out.println("g.all = " + g.all);
//      System.out.println("rest: " + g.all.subtract(vertices));
//      new Subgraph(g, g.all.subtract(vertices)).h.printRaw(System.out);
//      System.out.println("g has " + g.n + " vertices and " + 
//          g.numberOfEdges() + " edges");
      
      for (XBitSet compo: components) {
        XBitSet nb = g.neighborSet(compo);
        XBitSet clique = nb.convert(conv);
        h.fill(clique);
        cliques.add(clique);
        cliqueMap.put(clique, compo);
      }
    }
  }




  public void printVertices() {
    System.out.println("LocalGraph Vertices and Mapping:");
    System.out.println("Original Graph Vertices (inv):");
    for (int i = 0; i < inv.length; i++) {
      System.out.println("  Local Vertex " + i + " -> Original Vertex " + inv[i]);
    }
  }
}
