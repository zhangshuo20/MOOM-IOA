$title Epsilon Constraint Method for mop Model

$onText
The eps-Constraint Method

Additional information can be found at:

http://www.gams.com/modlib/adddocs/epscm.pdf

Mavrotas, G, Effective implementation of the eps-constraint method in
Multi-Objective Mathematical Programming problems.
Applied Mathematics and Computation 213, 2 (2009), 455-465.

Keywords: linear programming, eps-constraint method, multiobjective optimization
$offText

$if %system.filesys% == UNIX $abort.noerror 'This model cannot run on a non-Windows platform';
$call msappavail -Excel
$if errorlevel 1 $abort.noerror 'Microsoft Excel is not available!';

$inlineCom [ ]
$eolCom //

$stitle  Model Definitions

Sets
   i     'Index of sector' / N_1*N_48 /
   t     'Index of time period'/t_1*t_11/
   ie(i) 'index of electricity sectors' /N_25*N_30/
   k     'objective functions'   / economy, environment, employment, energy/;

Alias (i,j);

$set min -1
$set max +1

Parameter dir(k) 'direction of the objective functions'
                 / economy %max%, environment %min%, employment %max%, energy %min% /;

Scalar
   rSector1   'the minimum of sectoral growth rate'   /0.8/
   rSector2   'the maxmum of sectoral growth rate'    /1.2/;

*Table      A_2017(i,j)     'input-output coefficients'
*Table      elec_2017(i,j)  'input-output coefficients in electricity sectors'

Parameters
   A_2017(i,j)        'input-output coefficients'
   elec_2017(j,i)     'input-output coefficients in electricity sectors'
   F_2020(j)          'final demands of sector i in 2017'
   CI(t,j)            'carbon emissions intensity of sector i in year t'
   EI(t,j)            'energy intensity of sector i in year t'
   LI(t,j)            'employment intensity of sector i'
   output_2020(j)     'sectoral output of sector i in 2017'
   elec_demand(t)     'electricity demand in year t'
   lim_energy(t)      'limit of energy consumption in year t'
   lim_emp(t)
   iwater(t,j)         'water use intensity'
   iSO2(t,j)
   iNOX(t,j)
   iSD(t,j)
   iCOD(t,j)
   iAN(t,j)
   rGDP (t)           'the minimum of GDP growth rate in year t'
   x0(j)              'base (2017) year output of sector j'
   x2017(j)
   inco(j)            'intermediate coefficient'
   einco(j)           'intermediate coefficient of electricity sectors'
   opt(t)             'total output in each year'
   env(t)
   emp(t)
   ene(t)
   d_elec(t)

;
$CALL GDXXRW tst.xlsx Index=Index!A1:E21 trace=0
$GDXIN tst.gdx
$LOAD A_2017 = A_2017
$LOAD F_2020 = F_2020
$LOAD elec_2017 = elec_2017
$LOAD CI = CI
$LOAD EI = EI
$LOAD LI = LI
$LOAD output_2020 = output_2020
$LOAD elec_demand = elec_demand
$LOAD lim_energy = lim_energy
$LOAD lim_emp = lim_emp
$LOAD iwater = iwater
$LOAD iSO2 = iSO2
$LOAD iNOX = iNOX
$LOAD iSD = iSD
$LOAD iCOD = iCOD
$LOAD iAN = iAN
$LOAD rGDP = rGDP
$LOAD x0 = x0
$LOAD x2017 = x2017
$GDXIN
;

inco(j) = sum(i,A_2017(i,j));
einco(j)= sum(i,elec_2017(j,i));

*display inco,einco;

Variable
     z(k)         'objective function variables';

Variable

     z1      'economy objective function variable'
     z2      'environment objective function variable'
     z3      'employment objective function variable'
     z4      'energy objective function variable'
;

Positive variable
     x(t,j)  'output of sectors' ;


*Initial value:
x.lo(t,j) = x0(j);

Equation
* constraints
    IO_balance1(i)
    IO_balance2(i)
    IO_balance3(i)
    IO_balance4(i)
    IO_balance5(i)
    IO_balance6(i)
    IO_balance7(i)
    IO_balance8(i)
    IO_balance9(i)
    IO_balance10(i)
    IO_balance11(i)
    x2020(j)
    xlow(t,j)
    xhigh(t,j)
    GDP_growth(t)
    d_energy(t)
    t_emp(t)
* objective function
    objeco   'objective for maximizing economic outputs'
    objenv   'objective for minimizing carbon emissions'
    objemp   'objective for maximizing employment'
    objene   'objective for minimizing energy'
;

* constraints
IO_balance1(i) ..  x('t_1',i)-sum(j,A_2017(i,j)*x('t_1',j))=g=F_2020(i);
IO_balance2(i) ..  x('t_2',i)-sum(j,A_2017(i,j)*x('t_2',j))=g=F_2020(i);
IO_balance3(i) ..  x('t_3',i)-sum(j,A_2017(i,j)*x('t_3',j))=g=F_2020(i);
IO_balance4(i) ..  x('t_4',i)-sum(j,A_2017(i,j)*x('t_4',j))=g=F_2020(i);
IO_balance5(i) ..  x('t_5',i)-sum(j,A_2017(i,j)*x('t_5',j))=g=F_2020(i);
IO_balance6(i) ..  x('t_6',i)-sum(j,A_2017(i,j)*x('t_6',j))=g=F_2020(i);
IO_balance7(i) ..  x('t_7',i)-sum(j,A_2017(i,j)*x('t_7',j))=g=F_2020(i);
IO_balance8(i) ..  x('t_8',i)-sum(j,A_2017(i,j)*x('t_8',j))=g=F_2020(i);
IO_balance9(i) ..  x('t_9',i)-sum(j,A_2017(i,j)*x('t_9',j))=g=F_2020(i);
IO_balance10(i)..  x('t_10',i)-sum(j,A_2017(i,j)*x('t_10',j))=g=F_2020(i);
IO_balance11(i)..  x('t_11',i)-sum(j,A_2017(i,j)*x('t_11',j))=g=F_2020(i);

x2020(j) ..       x('t_1',j)=e=output_2020(j);

xlow(t-1,j) ..     x(t,j)=g=rSector1*x(t-1,j);

xhigh(t-1,j) ..    x(t,j)=l=rSector2*x(t-1,j);

GDP_growth(t-1) .. sum(j,x(t,j)-x(t,j)*inco(j))=g=(1+rGDP(t))*sum(j,x(t-1,j)-x(t-1,j)*inco(j));

d_energy(t) ..     sum(j,EI(t,j)*(x(t,j)-x(t,j)*inco(j)))=l=lim_energy(t);

t_emp(t) ..        sum(j,LI(t,j)*(x(t,j)-x(t,j)*inco(j)))=l=lim_emp(t);


* objective function

objeco ..   z('economy') =e= sum((j,t),x(t,j)-x(t,j)*inco(j));
objenv ..   z('environment') =e= sum((j,t),CI(t,j)*(x(t,j)-x(t,j)*inco(j)));
objemp ..   z('employment') =e=  sum((j,t),LI(t,j)*(x(t,j)-x(t,j)*inco(j)));
objene ..   z('energy') =e= sum((j,t),EI(t,j)*(x(t,j)-x(t,j)*inco(j)));

Model multiobj / all /;
option optCR = 0;


$STitle eps-Constraint Method
Set
   k1(k)  'the first element of k'
   km1(k) 'all but the first elements of k';

k1(k)$(ord(k) = 1) = yes;
km1(k)  = yes;
km1(k1) =  no;

Set kk(k) 'active objective function in constraint allobj';

Parameter
   rhs(k)     'right hand side of the constrained obj functions in eps-constraint'
   maxobj(k)  'maximum value from the payoff table'
   minobj(k)  'minimum value from the payoff table';

Variable
   a_objval   'auxiliary variable for the objective function'
   obj        'auxiliary variable during the construction of the payoff table';

Positive Variable
   sl(k)      'slack or surplus variables for the eps-constraints';

Equation
   con_obj(k) 'constrained objective functions'
   augm_obj   'augmented objective function to avoid weakly efficient solutions'
   allobj     'all the objective functions in one expression';

con_obj(km1).. z(km1) - dir(km1)*sl(km1) =e= rhs(km1);

* We optimize the first objective function and put the others as constraints
* the second term is for avoiding weakly efficient points
augm_obj.. sum(k1,dir(k1)*z(k1)) + 1e-3*sum(km1,sl(km1)/(maxobj(km1) - minobj(km1))) =e= a_objval;

allobj..   sum(kk, dir(kk)*z(kk)) =e= obj;

Model
   mod_payoff    / multiobj, allobj /
   mod_epsmethod / multiobj, con_obj, augm_obj /;

option limrow=0, limcol=0;
option solprint=off, solvelink=%solvelink.CallModule%;


option limRow = 0, limCol = 0, solPrint = off, solveLink = %solveLink.CallModule%;

Parameter payoff(k,k) 'payoff tables entries';

Alias (k,kp);

* Generate payoff table applying lexicographic optimization
loop(kp,
   kk(kp) = yes;
   repeat
      solve mod_payoff using LP maximizing obj;
      payoff(kp,kk) = z.l(kk);
      z.fx(kk)      = z.l(kk); // freeze the value of the last objective optimized
      kk(k++1)      = kk(k);   // cycle through the objective functions
   until kk(kp);
   kk(kp) = no;
*  release the fixed values of the objective functions for the new iteration
   z.up(k) =  inf;
   z.lo(k) = -inf;
);

if(mod_payoff.modelStat <> %modelStat.optimal% and
   mod_payoff.modelStat <> %modelStat.feasibleSolution%,
   abort 'no feasible solution for mod_payoff');

display payoff;

minobj(k) = smin(kp,payoff(kp,k));
maxobj(k) = smax(kp,payoff(kp,k));

$set fname p.%gams.scrext%

File fx 'solution points from eps-method' / "%gams.scrdir%%fname%" /;

$if not set gridpoints $set gridpoints 10

Set
   g         'grid points' / g0*g%gridpoints% /
   grid(k,g) 'grid';

Parameter
   gridrhs(k,g) 'rhs of eps-constraint at grid point'
   maxg(k)      'maximum point in grid for objective'
   posg(k)      'grid position of objective'
   firstOffMax  'counter'
   lastZero     'counter'
   numk(k)      'ordinal value of k starting with 1'
   numg(g)      'ordinal value of g starting with 0';

lastZero = 1;
loop(km1,
   numk(km1) = lastZero;
   lastZero  = lastZero + 1;
);
numg(g) = ord(g) - 1;

grid(km1,g) = yes; // Here we could define different grid intervals for different objectives
maxg(km1)   = smax(grid(km1,g), numg(g));
gridrhs(grid(km1,g))$(%min% = dir(km1)) = maxobj(km1) - numg(g)/maxg(km1)*(maxobj(km1) - minobj(km1));
gridrhs(grid(km1,g))$(%max% = dir(km1)) = minobj(km1) + numg(g)/maxg(km1)*(maxobj(km1) - minobj(km1));
display gridrhs;

* Walk the grid points and take shortcuts if the model becomes infeasible
posg(km1) = 0;
scalar  iter  'total number of iterations';
iter = 0;
*option savepoint=2;
repeat
   iter = iter + 1;
   rhs(km1) = sum(grid(km1,g)$(numg(g) = posg(km1)), gridrhs(km1,g));
   solve mod_epsmethod maximizing a_objval using LP;

   if(mod_epsmethod.modelStat <> %modelStat.optimal%,  // not optimal is in this case infeasible
      lastZero = 0;
      loop(km1$(posg(km1)  > 0 and lastZero = 0), lastZero = numk(km1));
      posg(km1)$(numk(km1) <= lastZero) = maxg(km1); // skip all solves for more demanding values of rhs(km1)
   else
      loop(k, put fx z.l(k):12:2); put /;
      Put_utility 'gdxout' / '..\results\main_eps2\result_'iter:0:0'.gdx';
      Execute_unload;
      Put_utility 'exec' / 'gdx2xls  ..\results\main_eps2\result_'iter:0:0'.gdx'
   );

*  Proceed forward in the grid
   firstOffMax = 0;
   loop(km1$(posg(km1) < maxg(km1) and firstOffMax = 0),
      posg(km1)   = posg(km1) + 1;
      firstOffMax = numk(km1);
   );
   posg(km1)$(numk(km1) < firstOffMax) = 0;
until sum(km1$(posg(km1) = maxg(km1)),1) = card(km1) and firstOffMax = 0;
putClose fx; // close the point file


* Get unique solutions from the point file using some Posix Tools (awk, (g)sort, uniq) that come with GAMS
$set awkscript awk.%gams.scrext%
file fa / "%gams.scrdir%%awkscript%" /; put fa 'BEGIN { printf("Table solutions(*,*)\n$ondelim\nsol';
loop(k, put ',' k.tl:0); putclose '\n"); }' / '{ print NR,$0 }' / 'END { print ";" }';
$if     %system.filesys% == UNIX execute 'cd "%gams.scrdir%" && sort %fname% | uniq | awk -f %awkscript% > g.%gams.scrext% && gams g.%gams.scrext% o=gx.%gams.scrext% lo=0 gdx=soleps';
$if NOT %system.filesys% == UNIX execute 'cd "%gams.scrdir%" && gsort %fname% | uniq | awk -f %awkscript% > g.%gams.scrext% && gams g.%gams.scrext% o=gx.%gams.scrext% lo=0 gdx=soleps';
execute 'mv -f "%gams.scrdir%soleps_e.gdx" .';

Set b Solutions /1*100/; Parameter solutions(b,k) Unique solutions;
execute_load 'soleps_e', solutions; display solutions;
$exit


