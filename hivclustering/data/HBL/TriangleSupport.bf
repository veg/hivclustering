TRY_NUMERIC_SEQUENCE_MATCH = 1;

first_call = 1;

function _testNetworkTriangle (filter, efv, seq1, seq2, seq3) {
// seq_i is the list of sequence indices within the datafilter
        
    if (first_call) {    
        global R1 = 2;
        global R2 = 2;
        qTN93 = {{*,t,R1*t,t}
                 {t,*,t,R2*t}
                 {R1*t,t,*,t}
                 {t,R2*t,t,*}};
             
        Model TN93 = (qTN93, efv, 1);
        Tree          T3 = (1,2,3);
        first_call = 0; 
     } else {
        R1 = 2;
        R2 = 2;
    }
    
    filterString = Join (",", {{seq1,seq2,seq3}});
    DataSetFilter D3 = CreateFilter (^filter, 1,, filterString);
    
    USE_LAST_RESULTS       = 0;
    OPTIMIZATION_METHOD    = 4;
        
    ClearConstraints (T3);
    LikelihoodFunction L3 = (D3,T3);
    Optimize (full, L3);
    
    
    // stash all MLEs
    
    GetString (_lfInfo,L3,-1);
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
    
    USE_LAST_RESULTS = 1;
    OPTIMIZATION_METHOD = 0;
    
    p_values = {3,1};
    bl = BranchLength (T3, -1);
    
    for (b = 0; b < 3; b+=1) {
        // restore MLEs 
        if (bl[b] > 0) {
            for (k = 0; k < Abs (MLE_stash); k+=1) {
                id = MLE_stash["INDEXORDER"][k];
                Eval ("`id` = " + MLE_stash[id]);
            }
            // set branch constraint
            ExecuteCommands ("ReplicateConstraint (\"this1.?:=0\",T3." + (b+1) + ");");
            Optimize (constrained, L3);
            p_values[b] = 0.5*(1-CChi2 (2*(full[1][0]-constrained[1][0]), 1));
        } else {
            p_values [b] = 0.5;
        }
    }
    
    return p_values;
           
}



map = {triangle_count, 3};

//fprintf (stdout, "Loading...\n");
DataSet       ds           = ReadDataFile (_py_sequence_file);
DataSetFilter filteredData = CreateFilter (ds,1);

nameToIndex = {};

//fprintf (stdout, "Indexing...\n");
for (k = 0; k < filteredData.species; k+=1) {
    GetString (seq_name, filteredData, k);
    nameToIndex[seq_name] = k + 1;
}

//fprintf (stdout, "Getting frequencies...\n");

COUNT_GAPS_IN_FREQUENCIES = 0;
HarvestFrequencies          (globalFreqs, filteredData, 1,1,1);


triangle_count = Abs(_py_triangle_sequences) $ 3;
all_p_values = {triangle_count, 3};


for (_t = 0; _t < triangle_count; _t += 1) {
    _toffset = _t * 3;
    
    //fprintf (stdout, _t, "\n");
    
    _sidx1 = nameToIndex[_py_triangle_sequences[_toffset]] - 1;
    assert (_sidx1 >= 0, "Failed to map " + _py_triangle_sequences[_toffset]);
    _sidx2 = nameToIndex[_py_triangle_sequences[_toffset + 1]] - 1;
    assert (_sidx2 >= 0, "Failed to map " + _py_triangle_sequences[_toffset + 1]);
    _sidx3 = nameToIndex[_py_triangle_sequences[_toffset + 2]] - 1;
    assert (_sidx3 >= 0, "Failed to map " + _py_triangle_sequences[_toffset + 2]);
    
    
    lpv = _testNetworkTriangle ("filteredData", globalFreqs,  
                                                    _sidx1,
                                                    _sidx2,
                                                    _sidx3);
                                                     
    for (z = 0; z < 3; z+=1) {                                                 
        all_p_values [_t][z] = lpv[z];                                                 
    } 
    
}


function _THyPhyAskFor(key)
{
    key_n = 0 + key;
    if (key_n >= 0 && key_n < triangle_count) {
    	return all_p_values[key_n][-1];
    }
    

    
    return "_THyPhy_NOT_HANDLED_";
}
