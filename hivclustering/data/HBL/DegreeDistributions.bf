harmonicNumbers    = {};
genHarmonicNumbers = {};

function logfactorial (k)
{
	f = 0;
	for (; k>1; k = k-1)
	{
		f += Log(k);
	}
	return f;
}

//--------------------------------------------------
// Start densities
//--------------------------------------------------


function waringDensity (x, rho, phi) {
	return Log(rho-1) + LnGamma (phi+rho) -LnGamma(phi+1) + LnGamma(x+phi) - LnGamma (x+phi+rho);
}

function yuleDensity (x, a) {
	return Log (a-1) + LnGamma (x) + LnGamma (a) - LnGamma (x+a);
}

function negBinomDensity (x, p, r) {
	return LnGamma (r+x) - LnGamma(r) + r*Log(1-p) + x*Log(p) - logfactorial(x);
}

function paretoDensity (x, p) {
	hsn = harmonicNumbers[p];
	if (hsn == 0)
	{
		vec = {1,10000} ["(_MATRIX_ELEMENT_COLUMN_+1)^(-p)"];
		hsn = +vec;
		harmonicNumbers[p] = +hsn;
	}
	return x^(-p)/hsn;
}

function LogZipfMandelbrotDensity (x, q, s) {
	return -s*Log(x+q);
}

//--------------------------------------------------
// end densities
//--------------------------------------------------

function likeFuncZipfMandelbrot (data, a, b)
{
	logL = 0;	
	maxD = 0;
	data["sumupZipfMandelbrot"][""];
	maxD = {1,maxD}["Exp(LogZipfMandelbrotDensity(_MATRIX_ELEMENT_COLUMN_+1,a,b))"];
	return logL - Log(+maxD);
}

function likeFuncWaring (data, a, b)
{
	logL = 0;	
	data["sumupWaring"][""];
	return logL;
}

function likeFuncPareto (data, a)
{
	logL = 0;	
	data["sumupPareto"][""];
	return logL;
}

function likeFuncYule (data, rho)
{
	logL = 0;	
	data["sumupYule"][""];
	return logL;
}

function likeFuncNB (data, a, b)
{
	logL = 0;	
	data["sumupNB"][""];
	return logL;
}

function sumupWaring (key, value)
{
	logL += (0+value) * waringDensity (0 + key, a, b);
	return 0;
}

function sumupZipfMandelbrot (key, value)
{
    xval = 0 + key;
    maxD = Max (maxD, xval);
	logL += (0+value) * LogZipfMandelbrotDensity (xval, a, b);
	return 0;
}

function sumupYule (key, value)
{
	logL += (0+value) * yuleDensity (0 + key, rho);
	return 0;
}

function sumupPareto (key, value)
{
	logL += (0+value) * Log(paretoDensity (0 + key, a));
	return 0;
}

function sumupNB (key, value)
{
	logL += value * negBinomDensity (0 + key, a, b);
	return 0;
}

totalNon0 = 0;	
maxD	  = 0;
degreeCounts = {};

//allDegs = {"0": 43, "1": 6, "2": 6, "3": 1, "4": 3};

for (k = 0; k < Abs(allDegs); k += 1)
{
	if (allDegs[k] > 0)
	{
		degreeCounts[k+1] += allDegs[k];
		totalNon0 += allDegs[k];
		maxD	   = Max (maxD, k+1);
	}
}


fprintf (stdout, maxD, "\n", degreeCounts, "\n", totalNon0, "\n");
fprintf (stdout, "Degree,PDF,CDF\n");

sum = 0;
for (k = 1; k <= maxD; k+=1)
{
	if (degreeCounts[k]>0)
	{
		p = degreeCounts[k] / totalNon0;
		sum += p;
		
		fprintf (stdout, Join (",",{{k, p, sum}}), "\n");
	}
}

x2 :> -1;
x1 :> 1;
x1 :< 10000;
x2 :< 10000;

_degree_fit_results = {};

x1 = 3;
x2 = 0;



VERBOSITY_LEVEL = 0;


Optimize (res, likeFuncWaring(degreeCounts,x1,x2));
fprintf (stdout, "\n\nWaring:\nLog(L) = ", res[1][0], "\nBIC = ", -res[1][0]*2 + Log(totalNon0) * res[1][1], "\nrho = ", x1, "\np = ", (x1-2)/(x1+x2-1), "\n");

x1MLE = x1;
x2MLE = x2;
critLevel = 1.92073*2;
FindRoot (zU, res[1][0] - likeFuncWaring(degreeCounts,x1,x2) - critLevel,x1, x1MLE, 10000);
FindRoot (zL, res[1][0] - likeFuncWaring(degreeCounts,x1,x2) - critLevel,x1, 0, x1MLE);
x1 = x1MLE;
x2 = x2MLE;

_degree_fit_results ["Waring"] = {"logL" : res[1][0], "BIC": -res[1][0]*2 + Log(totalNon0) * res[1][1], "rho": x1, "p": (x1-2)/(x1+x2-1), "x1" : x1, "x2" : x2, "rho_ci": {{zL__,zU__}}};




x2 := 0;
x1 :> 1;
Optimize (res, likeFuncYule(degreeCounts,x1));
fprintf (stdout, "\n\nYule:\nLog(L) = ", res[1][0], "\nBIC = ", -res[1][0]*2 + Log(totalNon0) * res[1][1], "\nrho = ", x1, "\n");

x1MLE = x1;
FindRoot (zU, res[1][0] - likeFuncYule(degreeCounts,x1) - critLevel,x1, x1MLE, 10000);
FindRoot (zL, res[1][0] - likeFuncYule(degreeCounts,x1) - critLevel,x1, 0, x1MLE);
x1 = x1MLE;

_degree_fit_results ["Yule"] = {"logL" : res[1][0], "BIC": -res[1][0]*2 + Log(totalNon0) * res[1][1], "rho": x1, "rho_ci": {{zL__,zU__}}};

x1 :< 1;
x1 :> 0;
x1 = 0.5;
x2 = 0.5;
x2 :> 0.000001;

Optimize (res, likeFuncNB(degreeCounts,x1,x2));
fprintf (stdout, "\n\nNeg. Binomial:\nLog(L) = ", res[1][0], "\nBIC = ", -res[1][0]*2 + Log(totalNon0) * res[1][1], "\np = ", x1, "\nr = ", x2, "\n");
_degree_fit_results ["Negative Binomial"] = {"logL" : res[1][0], "BIC": -res[1][0]*2 + Log(totalNon0) * res[1][1], "p": x1};

x1 :> 0;
x1 :< 1e26;

Optimize (res, likeFuncPareto(degreeCounts,x1));
fprintf (stdout, "\n\nPareto:\nLog(L) = ", res[1][0], "\nBIC = ", -res[1][0]*2 + Log(totalNon0) * res[1][1], "\nrho = ", x1,"\n");
x1MLE = x1;
FindRoot (zU, res[1][0] - likeFuncPareto(degreeCounts,x1) - critLevel,x1, x1MLE, 10000);
FindRoot (zL, res[1][0] - likeFuncPareto(degreeCounts,x1) - critLevel,x1, 0, x1MLE);
x1 = x1MLE;

_degree_fit_results ["Pareto"] = {"logL" : res[1][0], "BIC": -res[1][0]*2 + Log(totalNon0) * res[1][1], "rho": x1, "rho_ci": {{zL__,zU__}}};


function _THyPhyAskFor(key)
{
    if (key == "BestDistro")
    {
    	bestBIC = 1e100;
    	bestDistro = "";
        for (k = 0; k < Abs (_degree_fit_results); k += 1) {
        	key = _degree_fit_results["INDEXORDER"][k];
        	if ((_degree_fit_results[key])["BIC"] < bestBIC) {
        		bestBIC = (_degree_fit_results[key])["BIC"];
        		bestDistro = key;
        	}
        }
        return bestDistro;
    }
    
    if ((key $ "_BIC$" )[0] > 0)
    {
        return (_degree_fit_results[key[0][Abs(key)-5]])["BIC"];
    }
    
    if (key == "Waring_p") {
        return (_degree_fit_results[key[0][Abs(key)-3]])["p"];
    }
    
    if (key == "Waring_PDF") {
        res = {Abs(allDegs),1};
        return res["Exp(waringDensity(_MATRIX_ELEMENT_ROW_+1,"+(_degree_fit_results ["Waring"])["x1"]+","+(_degree_fit_results ["Waring"])["x2"]+"))"];
   }
    
    if (key == "Yule_PDF") {
        res = {Abs(allDegs),1};
        return res["Exp(yuleDensity(_MATRIX_ELEMENT_ROW_+1,"+(_degree_fit_results ["Yule"])["rho"] + "))"];
    }

    if (key == "Parto_PDF") {
        res = {Abs(allDegs),1};
        return res["Exp(paretoDensity(_MATRIX_ELEMENT_ROW_+1,"+(_degree_fit_results ["Parto_PDF"])["rho"] + "))"];
    }

    if (key == "Waring_rho_ci") {
        return (_degree_fit_results ["Waring"])["rho_ci"];
    }
    
    if (key == "Yule_rho_ci") {
        return (_degree_fit_results ["Yule"])["rho_ci"];
    }

    if (key == "Pareto_rho_ci") {
        return (_degree_fit_results ["Pareto"])["rho_ci"];
    }

    if (key == "Waring" || key == "Yule" || key == "Pareto")
    {
    	return (_degree_fit_results[key])["rho"];
    }
    
    return "_THyPhy_NOT_HANDLED_";
}
