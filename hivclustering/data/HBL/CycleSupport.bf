TRY_NUMERIC_SEQUENCE_MATCH = 1;

first_call = 1;

function _testNetworkCycle (filter, efv, seq1, seq2, seq3, seq4) {
// seq_i is the list of sequence indices within the datafilter

    if (first_call) {
        global R1 = 2;
        global R2 = 2;
        qTN93 = {{*,t,R1*t,t}
                 {t,*,t,R2*t}
                 {R1*t,t,*,t}
                 {t,R2*t,t,*}};

        Model TN93 = (qTN93, efv, 1);
        Tree          T1234 = ((1,2)N12,3,4);
        Tree          T2341 = ((2,3)N23,4,1);
        //Tree          T3412 = ((3,4)34,1,2);
        // same as T1234
        Tree          T4123 = ((4,1)N14,2,3);

        first_call = 0;
     } else {
        R1 = 2;
        R2 = 2;
    }

    filterString = Join (",", {{seq1,seq2,seq3,seq4}});
    DataSetFilter D4_1 = CreateFilter (^filter, 1,, filterString);
    DataSetFilter D4_2 = CreateFilter (^filter, 1,, filterString);
    DataSetFilter D4_3 = CreateFilter (^filter, 1,, filterString);
    
    
    p_values = {4,1};
        // 1-2 edge, 2-3 edge, 3-4 edge, 4-1 edge

    ClearConstraints (T1234);
    ClearConstraints (T2341);
    ClearConstraints (T4123);
    
    // handle 1-2-3-4 and 3-4-1-2 cases first
    
    LikelihoodFunction L4_1234 = (D4_1,T1234);
    USE_LAST_RESULTS       = 0;
    OPTIMIZATION_METHOD    = 4;
    Optimize (full_1234, L4_1234);
    LikelihoodFunction L4_2341 = (D4_2,T2341);
    Optimize (full_2341, L4_2341);
    LikelihoodFunction L4_4123 = (D4_3,T4123);
    Optimize (full_4123, L4_4123);
    
    best_ll = full_1234[1][0];
    best_ll = Max (best_ll, full_2341[1][0]);
    best_ll = Max (best_ll, full_4123[1][0]);
    
    GetString (_lfInfo,L4_1234,-1);
    MLE_stash = {};

	_paramList = _lfInfo["Global Independent"];
    for (k = 0; k < Columns (_paramList); k+=1) {
        id = _paramList[k];
        MLE_stash [id] = Eval (id);
    }
	_paramList = _lfInfo["Local Independent"];
    for (k = 0; k < Columns (_paramList); k+=1) {
        id = _paramList[k];
        MLE_stash [id] = Eval (id);
    }

    bl = BranchLength (T1234, -1);
    // will be in the order 1,2,N12,3,4
    //fprintf (stdout, best_ll, full_1234[1][0],full_2341[1][0],full_4123[1][0], "\n");

    if (Abs (full_1234[1][0] - best_ll) < 1) {
        if (bl[1] > 0 || bl[3] > 0) {
            T1234.2.t := 0;
            T1234.3.t := 0;
            Optimize (null_1234, L4_1234);
            p_values[3] = 1-CChi2 (2*(full_1234[1][0]-null_1234[1][0]), 2);
            T1234.2.t = MLE_stash["T1234.2.t"];
            T1234.3.t = MLE_stash["T1234.3.t"];
        } else {
            p_values[3] = 1;
        }

        if (bl[0] > 0 || bl[4] > 0) {
            T1234.1.t := 0;
            T1234.3.t := 0;
            Optimize (null_1234, L4_1234);
            p_values[1] = 1-CChi2 (2*(full_1234[1][0]-null_1234[1][0]), 2);
            T1234.4.t = MLE_stash["T1234.4.t"];
            T1234.1.t = MLE_stash["T1234.1.t"];
        } else {
            p_values[1] = 1;
        }
    }


    bl = BranchLength (T2341, -1);
    // 2,3,N23,4,1
    //fprintf (stdout, bl, "\n");
    USE_LAST_RESULTS = 1;
    OPTIMIZATION_METHOD = 0;
    if (Abs (full_2341[1][0] - best_ll) < 1) {
        if (bl[1] > 0 || bl[3] > 0) {
            T2341.3.t := 0;
            T2341.4.t := 0;
            Optimize (null_2341, L4_2341);
            p_values[0] = 1-CChi2 (2*(full_2341[1][0]-null_2341[1][0]), 2);
        } else {
            p_values[0] = 1;
        }
    }

    bl = BranchLength (T4123, -1);
    // 4,1,N41,2,3
    //fprintf (stdout, bl, "\n");
    USE_LAST_RESULTS = 1;
    OPTIMIZATION_METHOD = 0;
    if (Abs (full_4123[1][0] - best_ll) < 1) {
        if (bl[1] > 0 || bl[3] > 0) {
            T4123.1.t := 0;
            T4123.2.t := 0;
            Optimize (null_4123, L4_4123);
            p_values[2] = 1-CChi2 (2*(full_4123[1][0]-null_4123[1][0]), 2);
        } else {
            p_values[2] = 1;
        }
    }

    //fprintf (stdout, p_values, "\n");
    return p_values;
    

}



map = {quad_count, 4};

//NORMALIZE_SEQUENCE_NAMES = FALSE;

/*_py_sequence_dump = ">N4
CCTCAAATCACTCTTTGGCAACGACCCATCGTCACAATAAGGATAGGCGGGCAACTAAGGGAAGCTCTATTAGACACAAGAGCAGATGATACAGTATTAGAAGACATGAATTTGCCAGGAAGATGGAAACCAAAGATGATAGGGGGAATTGGAGGTTTTGTCAAAGTAAGACAGTATGATCAGATACCCATAGAAATCTGTGGACATAAAGTTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACGTAATTGGAAGAAATTTGATGACTCAGCTTGGTTGCACTTTAAATTTTCCCATTAGTCCTATTGAAACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAAGTCAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTGCAGAACTGGAAAAGGAAGGAAAAATTTCAAAAATTGGACCTGAAAATCCATACAACACTCCAATATTTGCTATAAGGAAAAAAAACAGCACGAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGACTTCTGGGAAGTTCAATTAGGAATACCGCATCCCGCAGGGTTACCAAAGAACAGGTCAGTAACAGTACTGGATGTGGGTGATGCATATTTTTCAGTTCCCTTAGATAAAGACTTCAGGAAGTATACTGCATTTACTATACCTAGTATAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTTCCACAGGGATGGAAAGGATCACCAGCAATATTCCAAAGTAGCATGACAAAAATCTTAGAGCCTTTCAGAAAATAAAACCCAGATGTAGTTATCTACCAATACATGGATGATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAACAAAAGTAGAGGAACTGAGACAACATCTGTTGAAGTGGGGAATTTACACACCAGACAAAAAACATCAGAAAGAACCTCCATTCCTTTGGATGGGTTATGAACTCCATCCAGATAAATGGACAGTACAGCCAATAGTGCTGCCAGAAAAAGACAGCTGGACTGTCAATGACATACAGAAGTTAGTGGGAAAATTGAATTGGGCAAGTCAGATTTACCCAGGGATTAAAGTAAAACAATTATGTAAACTCCTTAGAGGAACCAAATCACTAACAGAAGTAATACCACTAACAAAAGAGGCAGAGCTAGAGCTGGCAGAAAACAGGGAGATTCTAAAACAACCAGTACATGGAGTGTATTATGACCCATCAAAAGACTTAATAGCAGAATTACAGAAGCAGGAGCAAGGC
>N3
CCTCAAATCACTCTTTGGCAACGACCCATCGTCACAATAAAGATAGGAGGGCAACTAAGGGAAGCTCTATTAGACACAAGAGCAGATGATACAGTATTAGAAGACATAAATTTGCCAGGAAGATGGAAACCAAAGATGATAGGGGGAATTGGAGGTTTTGTCAAAGTAAGACAGTATGATCAGATACCCATAGAAATCTGTGGACATAAAGTTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACGTAATTGGAAGAAATCTGATGACTCAGCTTGGTTGCACTTTAAATTTTCCCATTAGTCCTATTGAAACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAGAGTCAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTGCAGAACTGGAAAAGGAAGGAAAAATTTCAAAAATTGGACCTGAAAATCCATACAACACTCCAATATTTGCTATAAGGAAAAAGAACAGCACGAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGACTTCTGGGAAGTTCAATTAGGAATACCGCATCCCGCAGGGTTACCAAAGAACAGGTCAGTAACAGTACTGGATGTGGGTGATGCATATTTTTCAGTTCCTTTAGATAAAGACTTCAGGAAGTACACTGCATTTACTATACCTAGTATAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTTCCACAGGGATGGAAAGGATCACCAGCAATATTCCAAAGTAGCATGACAAAAATCTTAGAGCCTTTCAGAAAATAAAACCCAGATATAGTTATCTACCAATACATGGATGATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAACAAAAGTAGAGGAACTGAGACAACATCTGTTGAAGTGGGGATTTTACACACCAGACAAAAAACATCAGAAAGAACCTCCATTCCTTTGGATGGGTTATGAACTCCATCCAGATAAATGGACAGTACAGCCAATAGTGCTGCCAGAAAAAGACAGCTGGACTGTCAATGACATACAGAAGTTAGTGGGAAAATTGAATTGGGCAAGTCAGATTTACCCAGGGATTAAAGTAAAACAATTATGTAAACTCCTTAGAGGAACCAAATCACTAACAGAAGTAATACCACTAACAAAAGAGGCAGAGCTAGAGCTGGCAGAAAACAGGGAGATTCTAAAACAACCAGTACATGGAGTGTATTATGACCCATCAAAAGACTTAATAGCAGAATTACAGAAGCAGGAGCAAGGC
>N2
CCTCAAATCACTCTTTGGCAACGACCCATCGTCACAATAAAGATAGGAGGGCAACTAAGGGAAGCTCTATTAGACACAGGAGCAGATGATACAGTATTAGAAGACATAAATTTGCCAGGAAGATGGAAACCAAAGATGATAGGGGGAATTGGAGGTTTTGTCAAAGTAAAACAGTATGATCAGATACCCATAGAAATCTGTGGACATAAAGTTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACGTAATTGGAAGAAATCTGATGACTCAGCTTGGTTGCACTTTAAATTTTCCCATTAGTCCTATTGAAACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAGAGTCAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTGCAGAACTGGAAAAGGAAGGAAAAATTTCAAAAATTGGACCTGAAAATCCATACAATACTCCAATATTTGCTATAAGGAAAAAGAACAGCACGAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGACTTCTGGGAAGTTCAATTAGGAATACCACATCCCGCAGGGTTACCAAAGAACAGGTCAGTAACAGTACTGGATGTGGGTGATGCATATTTTTCAGTTCCTTTAGATAAAGACTTCAGGAAGTACACTGCATTTACCATACCTAGTATAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTTCCACAGGGATGGAAAGGATCACCAGCAATATTCCAAAGTAGCATGACAAAAATCTTAGAGCCTTTCAGAAAATAAAACCCAGATATAGTTATCTACCAATACATGGATGATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAACAAAAGTAGAGGAACTGAGACAACATCTGTTGAAGTGGGGATTTTACACACCAGACAAAAAACATCAGAAAGAACCTCCATTCCTTTGGATGGGTTATGAACTCCATCCAGATAAATGGACAGTACAGCCAATAGTGCTGCCAGAAAAAGACAGCTGGACTGTCAATGACATACAGAAGTTAGTGGGAAAATTGAATTGGGCAAGTCAGATTTACCCAGGGATTAAAGTAAAGCAATTATGTAAACTCCTTAGAGGAACCAAATCACTAACAGAAGTAATACCACTAACAAAAGAGGCAGAGCTAGAGCTGGCAGAAAACAGGGAGATTCTAAAACAACCAGTACATGGAGTGTATTATGACCCATCAAAAGACTTAATAGCAGAATTACAGAAGCAGGAGCAAGGC
>N1
CCTCAAATCACTCTTTGGCAACGACCCATCGTCACAATAAAGATAGGAGGGCAACTAAGGGAAGCTCTATTAGACACAGGAGCAGATGATACAGTATTAGAAGACATAAGTTTGCCAGGAAGATGGAAACCAAAGATGATAGGGGGAATTGGAGGTTTTGTCAAAGTAAAACAGTATGATCAGATACCCATAGAAATCTGTGGACATAAAGTTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACGTAATTGGAAGAAATCTGATGACTCAGCTTGGTTGCACTTTAAATTTTCCCATTAGTCCTATTGAAACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAGAGTCAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTGCAGAACTGGAAAAGGAAGGAAAAATTTCAAAAATTGGACCTGAAAATCCATACAATACTCCAATATTTGCTATAAGGAAAAAGAACAGCACGAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGACTTCTGGGAAGTTCAATTAGGAATACCACATCCCGCAGGGTTACCAAAGAACAGGTCAGTAACAGTACTGGATGTGGGTGATGCATATTTTTCAGTTCCCTTAGATAAAGACTTCAGGAAGTACACTGCATTTACCATACCTAGTATAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTTCCACAGGGATGGAAAGGATCACCAGCAATATTCCAAAGTAGCATGACAAAAATCTTAGAGCCTTTCAGAAAACAAAACCCAGATATAGTTATCTACCAATACATGGATGATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAACAAAAGTAGAGGAACTGAGACAACATCTGTTGAAGTGGGGATTTTACACACCAGACAAAAAACATCAGAAAGAACCTCCATTCCTTTGGATGGGTTATGAACTCCATCCAGATAAATGGACAGTACAGCCAATAGTGCTGCCAGAAAAAGACAGCTGGACTGTCAATGACATACAGAAGTTAGTGGGAAAATTGAATTGGGCAAGTCAGATTTACCCAGGGATTAAAGTAAAGCAATTATGTAAACTCCTTAGAGGAACCAAATCACTAACAGAAGTAATACCACTAACAGAAGAGGCAGAGCTAGAGCTGGCAGAAAACAGGGAGATTCTAAAACAACCAGTACATGGAGTGTATTATGACCCATCAAAAGACTTAATAGCAGAATTACAGAAGCAGGAGCAAGGC

";*/

DataSet       ds           = ReadFromString (_py_sequence_dump);
DataSetFilter filteredData = CreateFilter (ds,1);

nameToIndex = {};

//fprintf (stdout, "Indexing...\n");
for (k = 0; k < filteredData.species; k+=1) {
    GetString (seq_name, filteredData, k);
    nameToIndex[seq_name] = k + 1;
}

//fprintf (stdout, nameToIndex, "\n");

//fprintf (stdout, "Getting frequencies...\n");

COUNT_GAPS_IN_FREQUENCIES = 0;
HarvestFrequencies          (globalFreqs, filteredData, 1,1,1);

//_py_triangle_sequences = {"0" : "222-3wc|05222003|pol|plasma|pool|1|UCSD|NULL|NULL", "1" : "222-447|11132003|pol|plasma|pool|1|ViroLogic|03-142073|NULL", "2" : "222-47d|07022004|pol|plasma|pool|1|ViroLogic|04_128904|NULL", "3" : "222-4d8|12052005|pol|plasma|pool|1|ViroLogic|05_160264|NULL"};
//_py_triangle_sequences = {"0" : "N1", "1" : "N2", "2" : "N3", "3" : "N4"};

quad_count = Abs(_py_triangle_sequences) $ 4;
all_p_values = {quad_count, 4};

for (_t = 0; _t < quad_count; _t += 1) {
    _toffset = _t * 4;

    //fprintf (stdout, _t, "\n");

    _sidx1 = nameToIndex[_py_triangle_sequences[_toffset]] - 1;
    assert (_sidx1 >= 0, "Failed to map " + _py_triangle_sequences[_toffset]);
    _sidx2 = nameToIndex[_py_triangle_sequences[_toffset + 1]] - 1;
    assert (_sidx2 >= 0, "Failed to map " + _py_triangle_sequences[_toffset + 1]);
    _sidx3 = nameToIndex[_py_triangle_sequences[_toffset + 2]] - 1;
    assert (_sidx3 >= 0, "Failed to map " + _py_triangle_sequences[_toffset + 2]);
    _sidx4 = nameToIndex[_py_triangle_sequences[_toffset + 3]] - 1;
    assert (_sidx3 >= 0, "Failed to map " + _py_triangle_sequences[_toffset + 3]);


    lpv = _testNetworkCycle ("filteredData", globalFreqs,
                                                    _sidx1,
                                                    _sidx2,
                                                    _sidx3,
                                                    _sidx4);

    for (z = 0; z < 4; z+=1) {
        all_p_values [_t][z] = lpv[z];
    }
}


function _THyPhyAskFor(key)
{
    key_n = 0 + key;
    if (key_n >= 0 && key_n < quad_count) {
    	return all_p_values[key_n][-1];
    }



    return "_THyPhy_NOT_HANDLED_";
}
