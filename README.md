Graph Isomorphism Experiments (Java) — README
===========================================

This repository contains a Java implementation for running graph isomorphism (GI)
experiments based on the methodology described in the paper:

  “Finding Graph Isomorphisms in Heated Spaces in Almost No Time”
  Sara Najem and Amer E. Mouawad
  https://arxiv.org/abs/2601.03787

This code is intended as a research/experimental tool (not a general-purpose GI library).

--------------------------------------------------------------------
IMPORTANT LIMITATIONS AND ASSUMPTIONS
--------------------------------------------------------------------

1) Graph size
-------------
Default settings are designed for CONNECTED graphs with at most 500 vertices.

- By default the program rejects inputs with > 500 vertices (configurable).
- Large graphs may lead to long runtimes / high memory usage because the method
  computes dense Laplacian eigendecompositions repeatedly.

2) Connectivity
---------------
All input graphs must be CONNECTED.

- If a graph is disconnected, the program terminates with failure.
- Preprocess disconnected datasets (e.g., choose a connected component) first.

3) Graph type
-------------
- Undirected, unweighted graphs only
- No self-loops
- No parallel edges

--------------------------------------------------------------------
REQUIREMENTS
--------------------------------------------------------------------

- Java 19 (the pom.xml targets source/target 19)
- Maven 3.8+ recommended

Dependencies are handled by Maven (JGraphT + EJML + ojalgo).

--------------------------------------------------------------------
PROJECT STRUCTURE
--------------------------------------------------------------------

Expected Maven layout:

  pom.xml
  src/main/java/iso/Main.java
  src/main/java/iso/GraphUtils.java
  src/main/java/iso/EigenUtils.java
  src/main/java/iso/SpectralUtils.java
  src/main/java/iso/CurvatureUtils.java
  src/main/java/iso/SignatureUtils.java
  src/main/java/iso/Constants.java

--------------------------------------------------------------------
BUILD FROM SOURCE (Maven)
--------------------------------------------------------------------

From the repo root:

  mvn clean package

This project uses the Maven Shade plugin to build a self-contained ("fat") JAR
with the entry point set to:

  iso.Main

After a successful build, look in:

  target/

Typical output JAR name (may vary by Maven version/config):

  target/iso-0.0.1-SNAPSHOT-shaded.jar

You can confirm by listing:

  ls target/*.jar

--------------------------------------------------------------------
RUNNING THE SHADED JAR
--------------------------------------------------------------------

Folder mode (batch processing)
------------------------------
The program walks a directory recursively, loads each graph file, generates a
randomly permuted copy, and attempts to certify isomorphism.

Example:

  java -jar target/iso-0.0.1-SNAPSHOT-shaded.jar \
    --GRAPH_FOLDER /path/to/graphs \
    --OUTPUT_LOG /path/to/out/OUTPUT_LOGS.csv \
    --SEED 12345

Notes:
- The tool writes 3 CSVs alongside OUTPUT_LOG by suffixing:
    *_SUCCESS.csv, *_FAIL.csv, *_EXCEPTION.csv
- Already-processed graphs are skipped by reading the first column of those CSVs.

Pair mode
---------
This Java implementation is currently oriented around folder processing via GRAPH_FOLDER.
If you want a strict "compare graph1 vs graph2" mode, you can add a small entry point that calls `solve()` for explicit paths.

--------------------------------------------------------------------
INPUT GRAPH FORMATS
--------------------------------------------------------------------

The reader auto-detects several common formats and ignores comment lines.

Supported comment prefixes:
  \#   %   c   C

1) DIMACS / bliss-style format (usually 1-indexed)
--------------------------------------------------
Example:

  c example graph
  p edge 4 3
  e 1 2
  e 2 3
  e 3 4

- Vertices are often 1..n; the loader applies a shift heuristic to convert to 0..n-1.

2) Edge list with explicit header (0-indexed)
---------------------------------------------
Example:

  5 4
  0 1
  1 2
  2 3
  3 4

- First line: n m
- Then exactly m edges

3) Flexible edge list (no header)
---------------------------------
Example:

  1 2
  2 3
  3 4

- The loader infers vertex count and attempts a 0/1-index shift heuristic.

--------------------------------------------------------------------
WHAT THE PROGRAM DOES (HIGH-LEVEL)
--------------------------------------------------------------------

For each input graph G:
1) Validates format and connectedness
2) Generates a permuted copy G' using a fixed seed
3) Optionally subdivides edges until spectral dimension estimation succeeds
4) Iteratively refines vertex equivalence classes using:
   - Laplacian spectrum
   - spectral dimension estimate
   - heat-kernel diagonal + fitted coefficients (quantized)
   - BFS-layered vertex signatures
5) If ambiguous classes remain, attaches “gadgets” (cliques/paths) to break symmetry
6) Stops when all original vertices are either uniquely identified (singletons)
   or belong to twin-classes that match across both graphs, then verifies the mapping.

--------------------------------------------------------------------
CONFIGURATION / CLI OVERRIDES
--------------------------------------------------------------------

This program supports a simple override system for static parameters.

Supported forms:
- --KEY=value
- --KEY value
- --no-KEY         (sets boolean KEY=false)
- --KEY            (sets boolean KEY=true)

Common parameters:
- --GRAPH_FOLDER PATH
- --OUTPUT_LOG PATH
- --SEED LONG
- --N_THREADS INT
- --COEFFICIENTS INT
- --MAX_ORIG_NODES INT
- --MAX_ORIG_EDGES INT
- --MAX_INPUT_BYTES LONG
- --MIN_VERTICES_AFTER_SUBDIVISION INT
- --MAX_VERTICES_AFTER_SUBDIVISION INT
- --ORIG_SUBD INT
- --ORIG_GADGET_SIZE INT
- --USE_TRIPLETS_ON_PAIR_FAILURE true|false
- --SKIP_TRIPLETS_ON_FIRST_ROUND true|false
- --DEBUG true|false

Numerical knobs:
- --EIGEN_DEFAULT_SCALE_OV DOUBLE
- --DS_EPS_OV DOUBLE
- --CURV_COEFF_DEFAULT_SCALE_OV DOUBLE

Example tuning command:

  java -jar target/iso-0.0.1-SNAPSHOT-shaded.jar \
    --GRAPH_FOLDER ./graphs \
    --OUTPUT_LOG ./out/OUTPUT_LOGS.csv \
    --MAX_ORIG_NODES 500 \
    --N_THREADS 8 \
    --COEFFICIENTS 6 \
    --SEED 12345 \
    --DEBUG true

--------------------------------------------------------------------
COMMON FAILURE REASONS
--------------------------------------------------------------------

- Graph is disconnected
- Graph exceeds MAX_ORIG_NODES or MAX_ORIG_EDGES
- Input format invalid (self-loops, malformed lines, duplicates depending on loader)
- Graph becomes too large after subdivision (MAX_VERTICES_AFTER_SUBDIVISION)
- Numerical issues (e.g., spectral dimension fitting fails for too-small graphs)

--------------------------------------------------------------------
DISCLAIMER
--------------------------------------------------------------------

This code is provided for research/experimental purposes.

No guarantees are made regarding:
- performance beyond the tested graph size range
- stability under modified parameters
