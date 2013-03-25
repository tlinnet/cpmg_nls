#pragma rtGlobals=1		// Use modern global access method.
#include  <Global Fit 2>

Function analyzeR1r()
	
	Variable Pvalue
	Variable i = 0
	String nextDataSet
	String dataSetList=""
	String dataFolder
	String plotName

	do
		nextDataSet = GetBrowserSelection(i)
		if (strlen(nextDataSet) == 0)
			break
		endif
		dataSetList += nextDataSet + ";"
		i += 1
	while(1)
	
	
	
	Make/O/D/N=(i) Presults
	Make/O/T/N=(i) dataName
	Make/O/D/N=(i) R1
	Make/O/D/N=(i) R2
	Make/O/D/N=(i) kEX
	Make/O/D/N=(i) phi

	i = 0
	do
		dataFolder = StringFromList(i,dataSetList,";")
		if (strlen(dataFolder) == 0)
			break
		endif
		plotName = GetLeafName(dataFolder)
		SetDataFolder dataFolder
		Wave R1r_coef,R1r_exch_coef
		Pvalue=FitAndCompare(tiltAngle,slockpwr,R1r_rates)
		dataName[i]=plotName
		Presults[i]=Pvalue
		if (Pvalue < 0.05)
			R1[i]=R1r_exch_coef[0]
			R2[i]=R1r_exch_coef[1]
			kEX[i]=R1r_exch_coef[2]
			phi[i]=R1r_exch_coef[3]
		else
			R1[i]=R1r_coef[0]
			R2[i]=R1r_coef[1]
		endif
		plotData(plotName)
		plot2(plotName,i)
		i += 1
	while(1)
	
	
		
End		

Function plotData(plotName)
	String plotName
	Display/K=1 R1r_rates vs tiltAngle
	AppendToGraph fitNoExch, fitExch1, fitExch2, fitExch3
	ModifyGraph tick=2,mirror=1,minor=1,standoff=0
	ModifyGraph axThick=0.5
	ModifyGraph width={Aspect,1.618}
	Label left "\\f02R\\B1\\F'Symbol'r\\F'Arial'\\M \\f00(s\\S-1\\M)"
	Label bottom "tilt angle (¡)"
	ModifyGraph lblMargin=12,lblLatPos=0
	SetAxis/A/N=1/E=1 left
	SetAxis bottom 0,180
	ModifyGraph mode(fitNoExch)=0,rgb(fitNoExch)=(0,0,0)
	ModifyGraph mode(fitExch1)=0,rgb(fitExch1)=(48896,49152,65280)
	ModifyGraph mode(fitExch2)=0,rgb(fitExch2)=(32768,32768,65280)
	ModifyGraph mode(fitExch3)=0,rgb(fitExch3)=(0,0,65280)
	ModifyGraph mode(R1r_rates)=3,marker(R1r_rates)=8
	ModifyGraph zColor(R1r_rates)={slockpwr,*,*,Rainbow,0}
	ErrorBars/T=0.5/L=0.5 R1r_rates Y,wave=(R1r_err,R1r_err)
	TextBox/C/N=text0/Z=1/F=0/A=MT/E=1/X=0/Y=0 plotName
	Print plotName
End

Function FitAndCompare(tiltAngle,omega1,R1r_rates)
//Fit R1r data to two functions (without and with an exchange contribution to R2)
//Selects which model to use by the F-test
	Wave tiltAngle,omega1,R1r_rates
	Variable ss1,df1,ss2,df2
	
	Wave Presults
	Variable model,P
	Make/N=2/D/O R1r_coef={2,23}
	Make/N=4/D/O R1r_exch_coef={1,40,10000,100000}
	Make/N=180/D/O fitNoExch,fitExch1,fitExch2,fitExch3
	
	FuncFit/X=1 R1r R1r_coef R1r_rates /X=tiltAngle /D
	ss1=V_chisq
	df1=V_npnts-(V_nterms-V_nheld)
	fitNoExch = R1r(R1r_coef,x)
	Duplicate/O W_sigma,errNoExch
		
	FuncFit/H="0000"/X=1 R1r_exch R1r_exch_coef R1r_rates /X={tiltAngle,omega1} /D
	ss2=V_chisq
	df2=V_npnts-(V_nterms-V_nheld)
	WaveStats/Q slockpwr
	fitExch1 = R1r_exch(R1r_exch_coef,x,1341.1)
	fitExch2 = R1r_exch(R1r_exch_coef,x,800.5)
	fitExch3 = R1r_exch(R1r_exch_coef,x,651.2)
	Duplicate/O W_sigma,errExch
	
	P=Ftest(ss1,df1,ss2,df2)
	return P
End

Function R1r_exch(w,x,y) : FitFunc
	Wave w
	Variable x
	Variable y

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x,y) = R1*cos(x*pi/180)^2+(R2+phi*kEX/((2*pi*y/tan(x*pi/180))^2+(2*pi*y)^2+kEX^2))*sin(x*pi/180)^2
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 2
	//CurveFitDialog/ x
	//CurveFitDialog/ y
	//CurveFitDialog/ Coefficients 4
	//CurveFitDialog/ w[0] = R1
	//CurveFitDialog/ w[1] = R2
	//CurveFitDialog/ w[2] = kEX
	//CurveFitDialog/ w[3] = phi

	return w[0]*cos(x*pi/180)^2+(w[1]+w[3]*w[2]/((2*pi*y/tan(x*pi/180))^2+(2*pi*y)^2+w[2]^2))*sin(x*pi/180)^2
End

Function Rex(w,x,y)
	Wave w
	Variable x,y
	
	return (w[3]*w[2]/((2*pi*y/tan(x*pi/180))^2+(2*pi*y)^2+w[2]^2))*sin(x*pi/180)^2
End

Function R1r(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = R1*cos(x*pi/180)^2+R2*sin(x*pi/180)^2
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 2
	//CurveFitDialog/ w[0] = R1
	//CurveFitDialog/ w[1] = R2

	return w[0]*cos(x*pi/180)^2+w[1]*sin(x*pi/180)^2
End

Function R2eff_and_plot()
	
	Variable i
	String nextDataSet
	String dataSetList=""
	String dataFolder
	String plotName
	Wave Presult
	
	
	do
		nextDataSet = GetBrowserSelection(i)
		if (strlen(nextDataSet) == 0)
			break
		endif
		dataSetList += nextDataSet + ";"
		i += 1
	while(1)
	
	i = 0
	
	do
		dataFolder = StringFromList(i,dataSetList,";")
		if (strlen(dataFolder) == 0)
			break
		endif
		plotName = GetLeafName(dataFolder)
		SetDataFolder dataFolder

		
		plot2(plotName,i)
		i += 1
		
	while(1)

End

Function plot2(plotName,i)
		String plotName
		Variable i
		Wave R1r_rates,tiltAngle,R1r_coef,R1r_exch_coef,R1r_err,errExch,errNoExch
		Wave Presults=root:Presults
		Duplicate/O R1r_rates,R2eff
		Duplicate/O R1r_rates,R2eff_err
		Make/O/D/N=200 R2eff_calc
		SetScale/I x 0,20000,"", R2eff_calc
	
		if (Presults[i] < 0.05)
			R2eff = R1r_rates/sin(tiltAngle/180*pi)^2- R1r_exch_coef[0]/tan(tiltAngle/180*pi)^2
			R2eff_err = sqrt((R1r_err/sin(tiltAngle/180*pi)^2)^2+(errExch[0]/tan(tiltAngle/180*pi)^2)^2)
			R2eff_calc = R1r_exch_coef[1]+R1r_exch_coef[3]*R1r_exch_coef[2]/((x*2*pi)^2+R1r_exch_coef[2]^2)	
		else
			R2eff = R1r_rates/sin(tiltAngle/180*pi)^2- R1r_coef[0]/tan(tiltAngle/180*pi)^2
			R2eff_err = sqrt((R1r_err/sin(tiltAngle/180*pi)^2)^2+(errNoExch[0]/tan(tiltAngle/180*pi)^2)^2)
			R2eff_calc = R1r_coef[1]	
		endif
		Display/K=1 R2eff vs w_eff
		ErrorBars R2eff Y,wave=(R2eff_err,R2eff_err)
		ModifyGraph mode=3
		ModifyGraph zColor(R2eff)={slockpwr,*,*,Rainbow,0}
		AppendToGraph R2eff_calc
		ModifyGraph marker(R2eff)=8
		TextBox/C/N=text0/Z=1/F=0/A=MT/E=1/X=0/Y=0 plotName
		//SetAxis bottom 0,10000
		SetAxis left 0,(R2eff_calc(0)+2)
		Label left "\\f02R\\f00\\B2,eff\\M (s\\S-1\\M)"
		Label bottom "\\F'Symbol'w\\B\\F'Arial'eff\\M (s\\S-1\\M)"
End
	


//	GetLeafName(pathStr)
//		Utility routine that returns the "leaf" of a path.
//		For example, if the path is "hd:Folder A:Folder B:File A", this returns "File A".
//		The path separator is assumed to be colon.
//		If no leaf name is found, returns "".
//		It is assumed that the input path does NOT end with a colon.
Function/S GetLeafName(pathStr)
	String pathStr
	
	String leafName
	Variable pos
	
	leafName = ""
	pos = strlen(pathStr) - 1								// Search from the end
	do
		if (CmpStr(pathStr[pos], ":") == 0)				// Found the colon?
			leafName = pathStr[pos+1, strlen(pathStr)]
			break											// Break out of loop.
		endif
		pos -= 1
	while(pos > 0)
	
	return leafName
End


//	LoadDataSetFromFile(pathName, fileName)
//		This function loads the waves from a sample data file into a new IGOR data folder.
//		It creates the data folder, giving it the name of the file.
//		It then loads the data from the file into the new data folder.
//
//		If pathName is "" or fileName is "", it will put up a dialog from which you can
//		choose the file to be loaded.
//		If both pathName and fileName are not "", then it will load the file without
//		putting up a dialog. In this case, pathName must contain the name of an IGOR
//		symbolic path that you have created which points to the folder containing the
//		sample data. fileName must be the name of a file within that folder.
Function LoadDataSetFromFile(pathName, fileName)
	String pathName			// Name of an existing IGOR symbolic path or ""
	String fileName				// Name of file within folder pointed to by symbolic path or ""

	// If either of the input parameters is "", put up a dialog to get file.
	if ((strlen(pathName)==0) %| (strlen(pathName)==0))
		Variable dummyRefNum
		String pathToFolder
		Open/R/D/M="Choose data file" dummyRefNum		// This sets an automatically created local variable named S_fileName.
		
		if (strlen(S_fileName) == 0)						// User cancelled?
			return -1
		endif
		
		// Now break the full path to file down into path to folder plus file name.
		// This is done by searching for the last colon in the full path.
		fileName = GetLeafName(S_fileName)
		if (strlen(fileName) == 0)
			Print "LoadDataSetFromFile bug"
			return -1
		endif
		pathToFolder = S_fileName[0, strlen(S_fileName) - strlen(fileName) - 1]
		
		// Make sure there is an IGOR symbolic path pointing to the folder
		NewPath/O/Q CurrentDataFilePath pathToFolder
		pathName = "CurrentDataFilePath"
	endif

	// Save the current data folder so we can restore it below.
	String savedDataFolder = GetDataFolder(1)
	
	// If it does not already exist, make a data folder to contain the waves from the file.
	// The new data folder is created in the current data folder with the name the
	// same as the name of the file we are about to load.
	// The /S flag sets new data folder as the current data folder.
	NewDataFolder/O/S :$fileName
	
	// This loads data set from the file into the current IGOR data folder
	// which is the data folder we just created above. The wave names come
	// from the first line of the data in the file because of the /W flag.
	LoadWave/G/D/O/A/Q/P=$pathName fileName
	Rename wave0,dof
	Rename wave1,slockpwr
	Rename wave2,w_eff
	Rename wave3,tiltAngle
	Rename wave4,R1r_rates
	Rename wave5,R1r_err

	
	// Restore the original data folder.
	SetDataFolder savedDataFolder
	
	return 0
End

//	LoadDataSetsFromFolder(pathName)
//		Loads data from all of the files in the folder associated with the specified symbolic path.
//		pathName is the name of an IGOR symbolic path which you can create with the New Path
//		dialog in the Misc menu.
//		If pathName is "", it puts up a dialog from which you can choose the folder.
//
//		NOTE: This function assumes that ALL of the files in the specified folder are data
//			    files and thus loads them all.
Function LoadDataSetsFromFolder(pathName)
	String pathName			// Name of an existing IGOR symbolic path or ""

	// If pathName is "", allow user to create a new symbolic path to the folder containing the runs of data.
	if (strlen(pathName) == 0)
		NewPath/O/Q/M="Choose folder containing data files" CurrentDataFilePath
		PathInfo CurrentDataFilePath		// Check to see if user created the path.
		if (V_flag == 0)
			return -1						// User cancelled.
		endif
		pathName = "CurrentDataFilePath"
	endif

	Variable i, numFilesLoaded
	String fileName
	
	numFilesLoaded = 0
	i = 0
	do
		fileName = IndexedFile($pathName, i, "????")	// This loads ALL files in the folder. See IndexedFile help for details.
		if (strlen(fileName) == 0)
			break
		endif
		LoadDataSetFromFile(pathName, fileName)
		i += 1
		numFilesLoaded += 1
	while(1)
	
	return numFilesLoaded
End

Function Ftest(ss1,df1,ss2,df2)
	Variable ss1,df1,ss2,df2
	Variable fRatio,Pvalue
	Variable a,b,x
	String conclusion
	
	fRatio = ((ss1-ss2)/(df1-df2))/(ss2/df2)
	a = df2/2
	b = (df1-df2)/2
	x = df2/(df2+(df1-df2)*fRatio)
	Pvalue = betai(a, b, x)
	return Pvalue
End
