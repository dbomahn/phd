using DelimitedFiles,DataFrames,JuMP,CPLEX,MathOptFormat,LinearAlgebra#,Statistics
# dir1 = "C:\\Users\\AK121396\\Desktop\\"
# dir2 = dir1*"varval\\PF\\"
dir1 = "/home/ak121396/Desktop/FLPInstances/FLPvlp"
dir2 = dir1*"/PF/"
# subdir = "/home/ak121396/Desktop/FLPInstances/size_5_10/"
# files = readdir(subdir)[4]
# f = readdlm(subdir*files)


#######################################################################
warming up cplex

cd("/home/ak121396/Downloads/KirlikSayin2014/")
run(`./main ./lpfiles/05_010_10.txt.lp`)
# calling files where var_vals are
# files = readdir(pwd()*"/ex/")
# f = readdlm(pwd()*"/ex/"*files[11], '\n' ,'\n',header=true)



##############KirlikSayin2014 ########

#include "Epsilon.h"


#include <fstream>
#include <iostream>
#include <string.h>
#include <numeric>

#include <sys/param.h>
#include <sys/times.h>
#include <sys/types.h>
#include <time.h>

using namespace std;

Epsilon::Epsilon(int _p) {
 counter = 0;
 number = 0;
 delta = 1;
 numeff = 0;

 //number of objective functions
 p = _p;

 //Cplex object
 cplex= IloCplex(env);

 //cplex parameters
 cplex.setParam(IloCplex::Param::MIP::Limits::Nodes, 1000000000);	//MIP node limit
 cplex.setParam(IloCplex::Param::MIP::Limits::TreeMemory, 1000000000);	//tree memory limit
 cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0.0);	//relative MIP gap tolerance
 cplex.setParam(IloCplex::Param::Simplex::Limits::Iterations, 1000000000); //absolute MIP gap iteration limit
 cplex.setParam(IloCplex::Param::ClockType, 1); 	//clock type: CPUtime
 cplex.setParam(IloCplex::Param::MIP::Display, 0);
 cplex.setParam(IloCplex::Param::Threads, 1); 		//number of threads
 cplex.setParam(IloCplex::Param::AdvInd, 1);
 //cplex.setParam(IloCplex::ParallelMode, -1);

 number = 0;
 delta = 1;

 //number of objective functions
 p = _p;

 //Set of objective functions
 objs.resize(p);
 for (i=0; i<p; i++) {
   objs[i]=IloExpr(env);
 }
 //Set of constraints
 cons=IloConstraintArray(env);

 //Resize variable array
 varArray.resize(p);
 //Resize num array
 numArray.resize(p);

 //all objective functions
 allobj = IloExpr(env);

 //lower and upper bound for the objective functions
 ideal = new int[p];
 upper = new int[p];

 //size of the T array
 sizet = (int) pow(2.0, p-1);

 //resize and allocate memory for the T_box
 T_box.resize(sizet);
 for (i=0; i<sizet; i++) {
   T_box[i] = new Box;
   T_box[i]->l = new int[p];
   T_box[i]->u = new int[p];
 }

}

Epsilon::~Epsilon(void) {
 //cout<<"destruct"<<endl;
 cplex.end();
 env.end();


 //deallocate memory T_box
 for (i=0; i<sizet; i++) {
   delete[] T_box[i]->l;
   delete[] T_box[i]->u;
   delete T_box[i];
 }
 T_box.clear();

 //deallocate
 delete[] ideal;
 delete[] upper;

}

void Epsilon::getParameters(char* fileName) {

 IloModel inputModel(env);
 IloExtractable extr;
 IloObjective objective;
 IloExpr expr;

 cplex.importModel(inputModel, fileName);

 j=0;
 for(IloModel::Iterator it(inputModel); it.ok(); ++it) {
   extr = it.operator *();

   //feasible region
   if (extr.isConstraint() == true) {
     //set of objective functions
     if (j<p) {
   IloRangeI* impl = dynamic_cast<IloRangeI *>((*it).asConstraint().getImpl());
   if (impl) {
     IloRange rng(impl);
     IloExpr expr = rng.getExpr();

     for (IloExpr::LinearIterator it2 = expr.getLinearIterator(); it2.ok(); ++it2) {
       objs[j] += it2.getCoef() * it2.getVar();
       varArray[j].push_back(it2.getVar());
       numArray[j].push_back(it2.getCoef());
     }
   }
   j=j+1;
     }
     //constraints
     else {
   cons.add(extr.asConstraint());
   number = number +1;
     }
   }
 }

 //all-objectives (for the second-stage)
 for (j=0; j<p; j++) {
   allobj += objs[j];
 }

 //create model object
 epsmodel = IloModel(env);

 //Initiliaze objective function
 objfunc = IloMinimize(env, allobj);
 epsmodel.add(objfunc);

 //Range for the objective functions
 range_obj.resize(p);
 for (j=0; j<p; j++) {
   range_obj[j] = IloRange(env, -IloInfinity, objs[j], IloInfinity);
     epsmodel.add(range_obj[j]);
 }

 //set of constraints (rows)
 for (i=0; i<number; i++) {
   epsmodel.add(cons[i]);
 }

 //extract the model
 cplex.extract(epsmodel);
}

int Epsilon::findMaxAreaBox() {
 //area is always positive
 float area = 0;
 //index of the largest area covered box
 int bindex = -1;

 for (i=0; i<(int)boxes.size(); i++) {
   if (boxes[i]->A != -1) {
     if (boxes[i]->A > area) {
   area = boxes[i]->A;
   bindex = i;
     }
   }
 }

 if (bindex == -1) {
   cout<< "Index = -1" <<endl;
   exit(0);
 }

 return bindex;

}


bool Epsilon::obtainBounds() {

 for (j=0; j<p; j++) {

   //set the objective function
   objfunc.setExpr(objs[j]);

   //solve the model. Since there is a feasible solution. this model always have solution
   cplex.solve();

   //check for the feasibility
   if (cplex.getCplexStatus() == CPX_STAT_INFEASIBLE )	{
     //model is infeasible
     return false;
   }

   //Ideal point for the j-th objective function
   ideal[j] = (int)floor(cplex.getObjValue() + 0.5);
 }


 //set the objective function as maximization
 objfunc.setSense(IloObjective::Maximize);

 for (j=0; j<p; j++) {
   //number of models solved
   counter = counter + 1;

   //set the objective function
   objfunc.setExpr(objs[j]);

   //solve the model. Since there is a feasible solution. this model always have solution
   cplex.solve();

   if (cplex.getCplexStatus() == CPXMIP_UNBOUNDED) {
     //cout<<"Max problem has an unbounded solution"<<endl;

     // set as infinity
     upper[j] = numeric_limits<int>::max();
   }
   else {
     //Obtain the upper point of the j-th objective function
     upper[j] = (int)floor(cplex.getObjValue() + 0.5);

     //there may be an efficient solution on the upper bounds of the j-th objective function
     upper[j] = upper[j] + 1;
   }
 }

 //set the objective function as minimization for the epsilon-constraint formulation
 objfunc.setSense(IloObjective::Minimize);

 //there exists a feasible solution to problem
 return true;

}

//objective function to minimize
bool Epsilon::epsilonModel(int k, const int* rhs, int* obj) {

 //increase the counter by 1.
 counter = counter + 1;

 //value of the variable
 double val;

 //Set the objective function as obj[k]
 objfunc.setExpr(objs[k]);

 for (j=0; j<p; j++) {
   if (j == k)	{
     range_obj[j].setUB(upper[j]);
   }
   else {
     range_obj[j].setUB(rhs[j] - delta);
   }
 }

 //solve the model
 cplex.solve();


 //solve model and check the feasibility of the model
 if (cplex.getCplexStatus() == IloCplex::Infeasible) {
   //model is infeasible
   return false;
 }
 else {
   //for the first lexicographic (get value of the k-th objective function)
   val = cplex.getObjValue() + 0.5;
   obj[k] = (int)floor(val);

   //add constraint for the k-th objective function
   range_obj[k].setUB(obj[k]);

   //first get objective values than add new constraint
   objfunc.setExpr(allobj);

   //solve the second stage
   cplex.solve();

   //get objective function values
   for (j=0; j<p; j++)	{
     if (j != k) {
   val = cplex.getValue(objs[j]) + 0.5;
   obj[j] = (int)floor(val);
     }
   }
   return true;
 }
}

void Epsilon::removeBox(int num, int* rhs, int* obj) {

 for (i=0; i<(int)boxes.size(); i++) {
   if (boxes[i]->A != -1) {
     found = true;
     for (j=0; j<p; j++) {
   if (j != num) {
     if ((boxes[i]->l[j] < obj[j]) || (boxes[i]->u[j] > rhs[j])) {
       found = false;
       break;
     }
   }
     }

     if (found == true) {
   // add to deleted boxes list
   liste.push_back(i);

   //delete the box
   boxes[i]->A = -1;
     }
   }
 }
}

bool Epsilon::checkEfficient(int* obj) {

 bool is_new = true;

 for (i=0; i < numeff; i++) {
   //check for all objective functions
   found = true;
   for (j=0; j<p; j++)	{
     if ( obj[j] - effset[i]->f[j] != 0) {
   found = false;
   break;
     }
   }

   if (found == true) {
     //sol is a previously obtained efficient solution
     is_new = false;

     break;
   }
 }

 if (is_new == true) {
   //create an object

   effset.resize(numeff + 1);
   effset[numeff] = new Efficient;
   effset[numeff]->f = new int[p];

   for (j=0; j<p; j++)	{
     effset[numeff]->f[j] = obj[j];
   }
   numeff = numeff + 1;

   return true;
 }
 else {
   //sol is a previously obtained efficient solution
   return false;
 }

}

void Epsilon::updateList(int num, int* obj) {
 int t, j1, tcount;

 //index list
 vector<int> tempL;
 //some of the boxes are appended to end of the list
 int lsize = (int)boxes.size();

 for (i=0; i<lsize; i++) {
   if (boxes[i]->A != -1) {
     //assign current box to T_box
     tcount = 1;
     for (j=0; j<p; j++) {
   if (j != num) {
     T_box[0]->l[j] = boxes[i]->l[j];
     T_box[0]->u[j] = boxes[i]->u[j];
   }
     }

     // Rectangular subdvision
     for (j=0; j<p; j++) {
   if (j != num) {
     if((boxes[i]->l[j] < obj[j]) &&  (obj[j] < boxes[i]->u[j])) {
       for (t=0; t < tcount; t++) {
         // subdivide the box into two new retangles
         for (j1=0; j1<p; j1++) {
       if (j1 != num) {
         if (j1 != j) {
           T_box[tcount + t]->l[j1] = T_box[t]->l[j1];
           T_box[tcount + t]->u[j1] = T_box[t]->u[j1];
         }
         else {
           T_box[tcount + t]->l[j] = obj[j];
           T_box[tcount + t]->u[j] = T_box[t]->u[j];
         }
       }
         }

         // All remain same except for the j-th upper vertex
         T_box[t]->u[j] = obj[j];
       }

       //update the tcount
       tcount = tcount * 2;
     }
   }
     }

     if (tcount > 1) {
   for (t=0; t < tcount; t++) {
     if (t == 0) {
       //initialize volume
       boxes[i]->A = 1.0;

       //assign last box to the current box (instead of removing box)
       for (j=0; j<p; j++)	{
         if (j != num) {
       boxes[i]->l[j] = T_box[t]->l[j];
       boxes[i]->u[j] = T_box[t]->u[j];
       boxes[i]->A = boxes[i]-> A * (T_box[t]->u[j] - ideal[j]);
         }
       }

       //take the logarithm of the rectangle volume to avoid huge numbers
       boxes[i]->A = log(boxes[i]->A);
     }
     else {
       if ((int)liste.size()>0) {
         sizeb = liste.front();
         liste.pop_front();

         tempL.push_back(sizeb);

         //initialize area
         boxes[sizeb]->A = -1;

         //assign last box to the current box (instead of removing box)
         for (j=0; j<p; j++) {
       if (j != num) {
         boxes[sizeb]->l[j] = T_box[t]->l[j];
         boxes[sizeb]->u[j] = T_box[t]->u[j];
       }
         }
       }
       else {
         sizeb = (int)boxes.size();
         boxes.resize(sizeb + 1);

         //allocate memmory
         boxes[sizeb] = new Box;
         boxes[sizeb]->l = new int[p];
         boxes[sizeb]->u = new int[p];
         boxes[sizeb]->A = 1.0;

         for (j=0; j<p; j++) {
       if (j != num) {
         boxes[sizeb]->l[j] = T_box[t]->l[j];
         boxes[sizeb]->u[j] = T_box[t]->u[j];
         boxes[sizeb]->A = boxes[sizeb]-> A * (T_box[t]->u[j] - ideal[j]);
       }
         }

         //take the logarithm of the rectangle volume
         boxes[sizeb]->A = log(boxes[sizeb]->A);
       }
     }
   }
     }
   }
 }

 for (i=0; i<(int)tempL.size(); i++) {
   //box list index
   sizeb = tempL[i];

   //initialize volume
   boxes[sizeb]->A = 1.0;

   for (j=0; j<p; j++) {
     if (j != num) {
   boxes[sizeb]->A = boxes[sizeb]-> A * (boxes[sizeb]->u[j] - ideal[j]);
     }
   }

   //take the logarithm of the rectangle volume
   boxes[sizeb]->A = log(boxes[sizeb]->A);
 }

}


void Epsilon::mainLoop(char* fileName) {

 //maximum volume rectangle index
 int bindex;

 //get parameters
 getParameters(fileName);

 bool check, feas;
 int* obj = new int[p];
 int* rhs = new int[p];

 //take the first objective function as an objective function of the epsilon-constraint model
 int num = 0;

 //obtain bounds for the model
 feas = obtainBounds();

 if (feas == false) {
   cout<<"Model has no feasible solution"<<endl;
   exit(0);
 }

 //extract model
 cplex.extract(epsmodel);

 //initial box
 boxes.resize(1);
 boxes[0] = new Box;
 boxes[0]->l = new int[p];
 boxes[0]->u = new int[p];
 boxes[0]->A = 1.0;

 //bounds for the problem
 for (j=0; j<p; j++) {
   if (j != num) {
     boxes[0]->l[j] = ideal[j];
     boxes[0]->u[j] = upper[j];
     boxes[0]->A = boxes[0]->A * (upper[j] - ideal[j]);
   }
 }

 //take the logarithm of the initial rectangle volume
 boxes[0]->A = log(boxes[0]->A);

 while (boxes.size() != liste.size()) {
   //get max sized area
   bindex = findMaxAreaBox();

   for (j=0; j<p; j++) {
     if (j != num) {
   rhs[j] = boxes[bindex]->u[j];
     }
   }

   //solve model with given rhs
   feas = epsilonModel(num, rhs, obj);

   //check for feasibility
   if (feas == true)	{
     //Solution is feasible

     //check the efficient solution whether obtain before
     check = checkEfficient(obj);

     if (check == true) {
   //a new efficient solution is obtained

   //update box list
   updateList(num, obj);

   //remove redundant boxes
   removeBox(num, rhs, obj);
     }
     else {
   //efficient solution has already been found
   removeBox(num, rhs, obj);
     }
   }
   else {
     //solution is infeasible
     removeBox(num, rhs, ideal);
   }
 }
 char* out = strcat(fileName,".ndf");
 ofstream effile(out, ofstream::out);
 for (i=0; i<(int)effset.size(); i++) {
   for (j=0; j<p; j++) {
     effile<<effset[i]->f[j]<<"\t";
   }
   effile<<endl;
 }
 effile.close();

 char* out2 = strcat(fileName,".log");
 ofstream logfile(out2, ofstream::out);
 logfile << counter;
 logfile.close();

 //clear memory

 //clear the list that stores all removed boxes
 liste.clear();

 //destruct boxes structure
 for (i=0; i<(int)boxes.size(); i++){
   delete[] boxes[i]->l;
   delete[] boxes[i]->u;
   delete boxes[i];
 }
 boxes.clear();

 // destruct nondominated set list
 for (i=0; i<(int)effset.size(); i++) {
   delete[] effset[i]->f;
   delete effset[i];
 }
 effset.clear();

 delete[] obj;
 delete[] rhs;
}
