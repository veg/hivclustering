LoadFunctionLibrary ("ReadDelimitedFiles");

DataSet ds = ReadDataFile (alignment);
GetString (seqNames, ds, -1);
cluster_definitions ["extractSequenceList"][""];

function extractSequenceList (key, value)
{
	sequenceIDs = {};
	splits = splitOnRegExp (value, ";");
	for (k = 0; k < Abs (splits); k+=1)
	{
		sequenceIDs [splits[k]] = 1;
	}
	
	DataSetFilter filteredData = CreateFilter (ds,1,"",sequenceIDs[seqNames[speciesIndex]]);
	DATA_FILE_PRINT_FORMAT = 4;
	writeTo = output_path + DIRECTORY_SEPARATOR + key + ".nex";
	//fprintf (stdout, writeTo, "\n");
	fprintf (writeTo, CLEAR_FILE, filteredData);
	return 0;
}
