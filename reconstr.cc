// This is inside a bigger file, but the rest has no interest at all.

G4double FunctionToReconstruction(const G4double *xx) {
	// Function I want to reconstruct, what is basically run the simulation and got the likelihood, called as "fun" (even if it is not funny because it does not work). 
	// I am sure that fun has a minimum, more or less clear.
	G4double fun;
	
	const G4double xMin = gposX;
	const G4double yMin = gposY;
	const G4double zMin = gposZ;
	const G4double tMin = xx[0];
	const G4double fMin = xx[1];

	PointSource(xMin,yMin,zMin,tMin,fMin);
	
	if (gsimevents > 0) {
		command.str("");
		command << "/run/beamOn " << gsimevents;
		UI->ApplyCommand(command.str());
	}
	
	fun = Minimization.Minimize();   // Minimization is a class I define with that name t build the likelihood function, not Minuit itself.
	G4cout <<  " !!!!!!!!!!!!!" <<fun << G4endl;
	return fun;
	
}

void ReconstructionOfEvents() {

	ROOT::Math::Minimizer* min = 
	      ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"); //This is the way I create the minimizer.
 
	min->SetMaxFunctionCalls(gReconstruction); //For some weird reason, this 3 parameters seem to don´t do anything. That is why I think Migrad fails...
	min->SetMaxIterations(10000);
	min->SetTolerance(0.001);
	
	ROOT::Math::Functor f(&FunctionToReconstruction,2); //Function to reconstruct with 2 dimensions
   
	// I define here 5 variables because at first I wanted to minimize a 5 dimensions function (lucky me). At the end I just use gtheta and gphi
	double variable[5] = {gposX,gposY,gposZ,gtheta,gphi};
	double step[5] = {10,10,10,1000,1000}; // This does nothing either, that is why now it has ridiculous values.
	double lower[5] = {200.,200.,200.,0,0};
	double upper[5] = {maxworld,maxworld,maxworld, 180, 360};
 
	min->SetFunction(f);

	min->SetLimitedVariable(0 , "tMin" , variable[3] , step[3] , lower[3], upper[3]);
	min->SetLimitedVariable(1 , "fMin" , variable[4] , step[4] , lower[4], upper[4]);

 
	min->Minimize(); 
	const double *xs = min->X();
	G4cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): " 
	    << FunctionToReconstruction(xs) << G4endl;  
 
	delete min;
}