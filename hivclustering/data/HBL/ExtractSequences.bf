
if (Type (_py_sequence_file) == "Matrix") {
    file_count = Rows(_py_sequence_file)*Columns (_py_sequence_file);

} else {
    file_count = Abs(_py_sequence_file);
}

DataSet       ds           = ReadDataFile (_py_sequence_file[0]);

name_to_id = {};

for (k = 1; k < file_count; k += 1) {
    DataSet read_more = ReadDataFile (_py_sequence_file[k]);
    DataSet ds = Combine (ds, read_more);
}

DataSetFilter filteredData = CreateFilter (ds,1);

for (k = 0; k < filteredData.species; k += 1) {
    GetString (seq_name, filteredData, k);
    name_to_id[seq_name] = k+1;
}


function _THyPhyAskFor(key) {

    if (key == "") {
        return filteredData.species;
    }

    if (name_to_id [key] > 0) {
        GetDataInfo (seqData, filteredData, name_to_id [key]-1);
        return seqData;
    }

    return "_THyPhy_NOT_HANDLED_";
}
