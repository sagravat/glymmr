# GLYMMR - GLYcan Motif MineR

This algorithm accepts a list of glycan structures in IUPAC format and outputs a set of substructures that represent frequently occurring "motifs". The input object has a **binder** and **nonbinder** element. Each element contains a list of **GlycanArrayData** objects that include the glycan id, linear code, and the RFU value. The output object is a list of objects that includes the motif, the reference of high binder **GlycanArrayData** objects that contain the motif, the reference of non binder **GlycanAraryData** objects that contain the motif, and an optional p-value and q-value using a One-Sample TTest corrected with a False Discovery Rate.

An alternative implementation using the Mann-Whitney Test has also been provided.