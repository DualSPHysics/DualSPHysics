/*===================================================================================
<MOORDYNPLUS> Copyright (c) 2025 by Ivan Martinez-Estevez

Ivan Martinez-Estevez and Jose M. Dominguez (Universidade de Vigo, Spain)
Matt Hall (github.com/mattEhall)

This file is part of MoorDynPlus. MoorDynPlus is free software: you can 
redistribute it and/or modify it under the terms of the GNU General Public 
License as published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

Linking the MoorDynPlus library statically or dynamically with other modules is
making a combined work based on this library. Thus, the terms and conditions
of the GNU General Public License cover the whole combination. As a special
exception, the copyright holders of MoorDynPlus give you permission to dynamically
link this library with the program DualSPHysics to produce a combined model
featuring the capabilities of both DualSPHysics and MoorDynPlus. This exception
is strictly limited to linking between the compiled MoorDynPlus library and
DualSPHysics. It does not extend to other programs or the use of the MoorDynPlus
source code beyond the stipulations of the GPL. When the exception is used,
this paragraph must be included in the copyright notice.

MoorDynPlus is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for details.

You should have received a copy of the GNU General Public License along with
MoorDynPlus. If not, see <http://www.gnu.org/licenses/>.
===================================================================================*/

/// \file IQSlines.cpp \brief Implements the class \ref IQSlines.

#include "IQSlines.h"

//====================================================================================================
// This function return the positions and tensions of a single mooring line
//====================================================================================================
int IQSlines::Catenary(double XF,double ZF,double L,double EA,double W,double CB,double Tol,
	double* HFout,double* VFout,double* HAout,double* VAout,int Nnodes,std::vector<double>& s,
	vector<double>& X,std::vector<double>& Z,std::vector<double>& Te,const unsigned lineid)
{																	// double s[10],double X[10],double Z[10],double Te[10])
	if(longwinded==1) cout << "In Catenary.  XF is " << XF << " and ZF is " << ZF << endl;

	double HF,VF,HA,VA;  // these are temporary and put to the pointers above with the "out" suffix

	// ===================== input/output variables =========================

	/*
	double XF; // in -Horizontal distance between anchor and fairlead (meters)
	double ZF; // in -Vertical   distance between anchor and fairlead (meters)
	double L;  // in -Unstretched length of line (meters)
	double EA; // in -Extensional stiffness of line (N)
	double W;  // in -Weight of line in fluid per unit length (N/m)
	double CB; // in -Coefficient of seabed static friction drag (a negative value indicates no seabed) (-)
	double Tol;// in -Convergence tolerance within Newton-Raphson iteration specified as a fraction of tension (-)
	double HF; // out-Effective horizontal tension in line at the fairlead (N)
	double VF; // out-Effective vertical tension in line at the fairlead (N)
	double HA; // out-Effective horizontal tension in line at the anchor   (N)
	double VA; // out-Effective vertical   tension in line at the anchor   (N)

	// can get rid of these node points???

	int N=10;  // in -Number of nodes where the line position and tension can be output (-)
	double s[N]; // in -Unstretched arc distance along line from anchor to each node where the line position and tension can be output (meters)
	double X[N]; // out-Horizontal locations of each line node relative to the anchor (meters)
	double Z[N]; // out-Vertical   locations of each line node relative to the anchor (meters)
	double Te[N];// out-Effective line tensions at each node (N)
	*/
	std::string warningText=""; ///<Stores the current warning to throw

	double LMax; // Maximum stretched length of the line with seabed interaction beyond which the line would have to double-back on itself; here the line forms an "L" between the anchor and fairlead (i.e. it is horizontal along the seabed from the anchor, then vertical to the fairlead) (meters)

	if(W>0.0) { // .TRUE. when the line will sink in fluid
		LMax=XF-EA/W+sqrt((EA/W)*(EA/W)+2.0*ZF*EA/W); // Compute the maximum stretched length of the line with seabed interaction beyond which the line would have to double-back on itself; here the line forms an "L" between the anchor and fairlead (i.e. it is horizontal along the seabed from the anchor, then vertical to the fairlead)
		if((L>=LMax) && (CB>=0.0)) {  // .TRUE. if the line is as long or longer than its maximum possible value with seabed interaction
			Log->PrintfWarning("Catenary: Unstretched line length too large %g m. The maximum stretched length is %g m. (line=%d)",L,LMax,lineid);
			if(longwinded==1) cout << "       d (horiz) is " << XF << " and h (vert) is " << ZF << " and L is " << L << endl;
			return -1;
		}
	}

	//- Initialize some commonly used terms that don't depend on the iteration:
	const double WL=W*L;
	const double WEA=W*EA;
	const double LOvrEA=L/EA;
	const double CBOvrEA=CB/EA;
	//- Smaller tolerances may take more iterations, so choose a maximum inversely proportional to the tolerance
	const int MaxIter=int(1.0/Tol);

	//- More initialization
	bool FirstIter=true; // Initialize iteration flag

	double dXFdHF; ///<Partial derivative of the calculated horizontal distance with respect to the horizontal fairlead tension (m/N): dXF(HF,VF)/dHF
	double dXFdVF; ///<Partial derivative of the calculated horizontal distance with respect to the vertical   fairlead tension (m/N): dXF(HF,VF)/dVF
	double dZFdHF; ///<Partial derivative of the calculated vertical   distance with respect to the horizontal fairlead tension (m/N): dZF(HF,VF)/dHF
	double dZFdVF; ///<Partial derivative of the calculated vertical   distance with respect to the vertical   fairlead tension (m/N): dZF(HF,VF)/dVF
	double EXF; ///<Error function between calculated and known horizontal distance (meters): XF(HF,VF)-XF
	double EZF; ///<Error function between calculated and known vertical   distance (meters): ZF(HF,VF)-ZF

	double DET; ///<Determinant of the Jacobian matrix (m^2/N^2)
	double dHF; ///<Increment in HF predicted by Newton-Raphson (N)
	double dVF; ///<Increment in VF predicted by Newton-Raphson (N)


	double HFOvrW; ///<= HF/W
	double HFOvrWEA; ///<= HF/WEA
	double Lamda0; ///<Catenary parameter used to generate the initial guesses of the horizontal and vertical tensions at the fairlead for the Newton-Raphson iteration (-)
	double LMinVFOvrW; ///<= L-VF/W
	double sOvrEA; ///<= s[I]/EA
	double SQRT1VFOvrHF2; ///<= SQRT( 1.0+VFOvrHF2      )
	double SQRT1VFMinWLOvrHF2; ///<= SQRT( 1.0+VFMinWLOvrHF2 )
	double SQRT1VFMinWLsOvrHF2; ///<= SQRT( 1.0+VFMinWLsOvrHF*VFMinWLsOvrHF )
	double VFMinWL; ///<= VF-WL
	double VFMinWLOvrHF; ///<= VFMinWL/HF
	double VFMinWLOvrHF2; ///<= VFMinWLOvrHF*VFMinWLOvrHF
	double VFMinWLs; ///<= VFMinWL+Ws
	double VFMinWLsOvrHF; ///<= VFMinWLs/HF
	double VFOvrHF; ///<= VF/HF
	double VFOvrHF2; ///<= VFOvrHF*VFOvrHF
	double VFOvrWEA; ///<= VF/WEA
	double Ws; ///<= W*s[I]
	double XF2; ///<= XF*XF
	double ZF2; ///<= ZF*ZF

				//- insertion-to get HF and VF initial guesses (FAST normally uses previous time step) 
	XF2=XF*XF;
	ZF2=ZF*ZF;

	if(L<=sqrt(XF2+ZF2)) { //.TRUE. if the current mooring line is taut
		Lamda0=0.2;
	}
	else { //The current mooring line must be slack and not vertical
		Lamda0=sqrt(3.0*((L*L-ZF2)/XF2-1.0));
	}
	// ! As above, set the lower limit of the guess value of HF to the tolerance
	HF=max(abs(0.5*W* XF/Lamda0),Tol);
	VF=0.5*W*(ZF/tanh(Lamda0)+L);

	//- end insertion ============

	/*
	! To avoid an ill-conditioned situation, ensure that the initial guess for
	!   HF is not less than or equal to zero.  Similarly, avoid the problems
	!   associated with having exactly vertical (so that HF is zero) or exactly
	!   horizontal (so that VF is zero) lines by setting the minimum values
	!   equal to the tolerance.  This prevents us from needing to implement
	!   the known limiting solutions for vertical or horizontal lines (and thus
	!   complicating this routine):
	*/

	HF=max(HF,Tol);
	XF=max(XF,Tol);
	ZF=max(ZF,Tol);

	// Solve the analytical, static equilibrium equations for a catenary (or
	//   taut) mooring line with seabed interaction:

	// Begin Newton-Raphson iteration:
	for(int I=1; I<=MaxIter; I++) {
		// Initialize some commonly used terms that depend on HF and VF:
		VFMinWL=VF-WL;
		LMinVFOvrW=L-VF/W;
		HFOvrW=HF/W;
		HFOvrWEA=HF/WEA;
		VFOvrWEA=VF/WEA;
		VFOvrHF=VF/HF;
		VFMinWLOvrHF=VFMinWL/HF;
		VFOvrHF2=VFOvrHF*VFOvrHF;
		VFMinWLOvrHF2=VFMinWLOvrHF*VFMinWLOvrHF;
		SQRT1VFOvrHF2=sqrt(1.0+VFOvrHF2);
		SQRT1VFMinWLOvrHF2=sqrt(1.0+VFMinWLOvrHF2);

		// Compute the error functions (to be zeroed) and the Jacobian matrix
		//   (these depend on the anticipated configuration of the mooring line):

		if((CB<0.0) || (W<0.0) || (VFMinWL>0.0)) { // .TRUE. when no portion of the line      rests on the seabed
			EXF=(log(VFOvrHF+SQRT1VFOvrHF2)-log(VFMinWLOvrHF+SQRT1VFMinWLOvrHF2))*HFOvrW+LOvrEA* HF-XF;
			EZF=(SQRT1VFOvrHF2-SQRT1VFMinWLOvrHF2)*HFOvrW+LOvrEA*(VF-0.5*WL)-ZF;
			dXFdHF=(log(VFOvrHF+SQRT1VFOvrHF2)-log(VFMinWLOvrHF+SQRT1VFMinWLOvrHF2))/W -
				((VFOvrHF+VFOvrHF2/SQRT1VFOvrHF2)/(VFOvrHF+SQRT1VFOvrHF2)
					- (VFMinWLOvrHF+VFMinWLOvrHF2/SQRT1VFMinWLOvrHF2)/(VFMinWLOvrHF+SQRT1VFMinWLOvrHF2))/W+LOvrEA;
			dXFdVF=((1.0+VFOvrHF/SQRT1VFOvrHF2)/(VFOvrHF+SQRT1VFOvrHF2)
				- (1.0+VFMinWLOvrHF/SQRT1VFMinWLOvrHF2)/(VFMinWLOvrHF+SQRT1VFMinWLOvrHF2))/W;
			dZFdHF=(SQRT1VFOvrHF2-SQRT1VFMinWLOvrHF2)/W-(VFOvrHF2/SQRT1VFOvrHF2-VFMinWLOvrHF2/SQRT1VFMinWLOvrHF2)/W;
			dZFdVF=(VFOvrHF/SQRT1VFOvrHF2-VFMinWLOvrHF/SQRT1VFMinWLOvrHF2)/W+LOvrEA;
		}
		else if(-CB*VFMinWL<HF) { // .TRUE. when a portion of the line rests on the seabed and the anchor tension is nonzero
			EXF=log(VFOvrHF+SQRT1VFOvrHF2)*HFOvrW-0.5*CBOvrEA*W* LMinVFOvrW*LMinVFOvrW+LOvrEA* HF+LMinVFOvrW-XF;
			EZF=(SQRT1VFOvrHF2-1.0)*HFOvrW+0.5*VF*VFOvrWEA-ZF;

			dXFdHF=log(VFOvrHF+SQRT1VFOvrHF2)/W-((VFOvrHF+VFOvrHF2/SQRT1VFOvrHF2)/(VFOvrHF+SQRT1VFOvrHF2))/W+LOvrEA;
			dXFdVF=((1.0+VFOvrHF/SQRT1VFOvrHF2)/(VFOvrHF+SQRT1VFOvrHF2))/W+CBOvrEA*LMinVFOvrW-1.0/W;
			dZFdHF=(SQRT1VFOvrHF2-1.0-VFOvrHF2/SQRT1VFOvrHF2)/W;
			dZFdVF=(VFOvrHF/SQRT1VFOvrHF2)/W+VFOvrWEA;
		}
		else {
			// 0.0<HF<=-CB*VFMinWL ! A portion of the line must rest on the seabed and the anchor tension is zero
			EXF=log(VFOvrHF+SQRT1VFOvrHF2)*HFOvrW-0.5*CBOvrEA*W*(LMinVFOvrW*LMinVFOvrW
				- (LMinVFOvrW-HFOvrW/CB)*(LMinVFOvrW-HFOvrW/CB))+LOvrEA* HF+LMinVFOvrW-XF;
			EZF=(SQRT1VFOvrHF2-1.0)*HFOvrW+0.5*VF*VFOvrWEA-ZF;
			dXFdHF=log(VFOvrHF+SQRT1VFOvrHF2)/W-((VFOvrHF+VFOvrHF2/SQRT1VFOvrHF2)/(VFOvrHF+SQRT1VFOvrHF2))/W
				+ LOvrEA-(LMinVFOvrW-HFOvrW/CB)/EA;
			dXFdVF=((1.0+VFOvrHF/SQRT1VFOvrHF2)/(VFOvrHF+SQRT1VFOvrHF2))/W+HFOvrWEA-1.0/W;
			dZFdHF=(SQRT1VFOvrHF2-1.0-VFOvrHF2/SQRT1VFOvrHF2)/W;
			dZFdVF=(VFOvrHF/SQRT1VFOvrHF2)/W+VFOvrWEA;
		}
		// Compute the determinant of the Jacobian matrix and the incremental
		//   tensions predicted by Newton-Raphson:

		DET=dXFdHF*dZFdVF-dXFdVF*dZFdHF;

		dHF=(-dZFdVF*EXF+dXFdVF*EZF)/DET;  //  ! This is the incremental change in horizontal tension at the fairlead as predicted by Newton-Raphson
		dVF=(dZFdHF*EXF-dXFdHF*EZF)/DET;  //  ! This is the incremental change in vertical   tension at the fairlead as predicted by Newton-Raphson

		dHF=dHF*(1.0-Tol*I);           // ! Reduce dHF by factor (between 1 at I=1 and 0 at I=MaxIter) that reduces linearly with iteration count to ensure that we converge on a solution even in the case were we obtain a nonconvergent cycle about the correct solution (this happens, for example, if we jump to quickly between a taut and slack catenary)
		dVF=dVF*(1.0-Tol*I);           // ! Reduce dHF by factor (between 1 at I=1 and 0 at I=MaxIter) that reduces linearly with iteration count to ensure that we converge on a solution even in the case were we obtain a nonconvergent cycle about the correct solution (this happens, for example, if we jump to quickly between a taut and slack catenary)

		dHF=max(dHF,(Tol-1.0)*HF);   //! To avoid an ill-conditioned situation, make sure HF does not go less than or equal to zero by having a lower limit of Tol*HF [NOTE: the value of dHF=( Tol-1.0 )*HF comes from: HF=HF+dHF=Tol*HF when dHF=( Tol-1.0 )*HF]
		// Check if we have converged on a solution, or restart the iteration, or
		//   Abort if we cannot find a solution:

		if((abs(dHF)<=abs(Tol*HF)) && (abs(dVF)<=abs(Tol*VF))) { // .TRUE. if we have converged; stop iterating! [The converge tolerance, Tol, is a fraction of tension]
			 //cout << "converged" << endl;
			break;
		}

		else if((I==MaxIter) && (FirstIter)) {  // .TRUE. if we've iterated MaxIter-times for the first time, try a new set of ICs;
			/*  ! Perhaps we failed to converge because our initial guess was too far off.
			!   (This could happen, for example, while linearizing a model via large
			!   perturbations in the DOFs.)  Instead, use starting values documented in:
			!   Peyrot, Alain H. and Goulois, A. M., "Analysis Of Cable Structures,"
			!   Computers & Structures, Vol. 10, 1979, pp. 805-813:
			! NOTE: We don't need to check if the current mooring line is exactly
			!       vertical (i.e., we don't need to check if XF==0.0), because XF is
			!       limited by the tolerance above.*/

			XF2=XF*XF;
			ZF2=ZF*ZF;

			if(L<=sqrt(XF2+ZF2)) { //.TRUE. if the current mooring line is taut
				Lamda0=0.2;
			}
			else { //The current mooring line must be slack and not vertical
				Lamda0=sqrt(3.0*((L*L-ZF2)/XF2-1.0));
			}

			HF=max(abs(0.5*W* XF/Lamda0),Tol); // ! As above, set the lower limit of the guess value of HF to the tolerance
			VF=0.5*W*(ZF/tanh(Lamda0)+L);

			// Restart Newton-Raphson iteration:
			I=0;
			FirstIter=false;
			dHF=0.0;
			dVF=0.0;
		}

		else if((I==MaxIter) && (!FirstIter)) { // .TRUE. if we've iterated as much as we can take without finding a solution; Abort
			Log->PrintfWarning("Warning from Catenary:: Iteration not convergent. Cannot solve quasi-static mooring line solution. Check ILine %d.",lineid);
			return -1;
		}

		// Increment fairlead tensions and iterate...
		HF=HF+dHF;
		VF=VF+dVF;
	}


	/* ! We have found a solution for the tensions at the fairlead!

	! Now compute the tensions at the anchor and the line position and tension
	!   at each node (again, these depend on the configuration of the mooring
	!   line):*/

	if((CB<0.0) || (W<0.0) || (VFMinWL>0.0)) { // .TRUE. when no portion of the line rests on the seabed
													  // Anchor tensions:
		HA=HF;
		VA=VFMinWL;
		//! ILine position and tension at each node:
		for(int I=0; I<Nnodes; I++) { // Loop through all nodes where the line position and tension are to be computed
			if((s[I]<0.0) || (s[I]>L)) {
				Log->PrintfWarning("Warning from Catenary:: All line nodes must be located between the anchor and fairlead (inclusive) in routine Catenary(). Check ILine %d.",lineid);
				//cout << "        s[I]=" << s[I] << " and L=" << L << endl;
				return -1;
			}
			Ws=W*s[I];    // Initialize
			VFMinWLs=VFMinWL+Ws;     // some commonly
			VFMinWLsOvrHF=VFMinWLs/HF;      // used terms
			sOvrEA=s[I]/EA;      // that depend on s[I]
			SQRT1VFMinWLsOvrHF2=sqrt(1.0+VFMinWLsOvrHF*VFMinWLsOvrHF);

			X[I]=(log(VFMinWLsOvrHF+SQRT1VFMinWLsOvrHF2)-log(VFMinWLOvrHF+SQRT1VFMinWLOvrHF2))*HFOvrW+sOvrEA* HF;
			Z[I]=(SQRT1VFMinWLsOvrHF2-SQRT1VFMinWLOvrHF2)*HFOvrW+sOvrEA*(VFMinWL+0.5*Ws);
			Te[I]=sqrt(HF*HF+VFMinWLs*VFMinWLs);
		}
	}
	else if(-CB*VFMinWL<HF) { // .TRUE. when a portion of the line rests on the seabed and the anchor tension is nonzero
								 // Anchor tensions:
		HA=HF+CB*VFMinWL;
		VA=0.0;

		// ILine position and tension at each node:
		for(int I=0; I<Nnodes; I++) {  // Loop through all nodes where the line position and tension are to be computed
			if((s[I]<0.0) || (s[I]>L)) {
				Log->PrintfWarning("Warning from Catenary:: All line nodes must be located between the anchor and fairlead (inclusive) in routine Catenary(). Check ILine %d.",lineid);
				//cout << "        s[I]=" << s[I] << " and L=" << L << endl;
				return -1;
			}

			Ws=W*s[I]; // ! Initialize
			VFMinWLs=VFMinWL+Ws; // ! some commonly
			VFMinWLsOvrHF=VFMinWLs/HF; // ! used terms
			sOvrEA=s[I]/EA; // ! that depend
			SQRT1VFMinWLsOvrHF2=sqrt(1.0+VFMinWLsOvrHF*VFMinWLsOvrHF); // ! on s[I]

			if(s[I]<=LMinVFOvrW) { // .TRUE. if this node rests on the seabed and the tension is nonzero
				X[I]=s[I]+sOvrEA*(HF+CB*VFMinWL+0.5*Ws*CB);
				Z[I]=0.0;
				Te[I]=HF+CB*VFMinWLs;
			}
			else {   // LMinVFOvrW<s<=L ! This node must be above the seabed
				X[I]=log(VFMinWLsOvrHF+SQRT1VFMinWLsOvrHF2)*HFOvrW+sOvrEA* HF+LMinVFOvrW-0.5*CB*VFMinWL*VFMinWL/WEA;
				Z[I]=(-1.0+SQRT1VFMinWLsOvrHF2)*HFOvrW+sOvrEA*(VFMinWL+0.5*Ws)+0.5* VFMinWL*VFMinWL/WEA;
				Te[I]=sqrt(HF*HF+VFMinWLs*VFMinWLs);
			}
		}
	}
	else { // 0.0< HF <=-CB*VFMinWL   !             A  portion of the line must rest  on the seabed and the anchor tension is    zero
		   // Anchor tensions:
		HA=0.0;
		VA=0.0;
		// ILine position and tension at each node:
		for(int I=0; I<Nnodes; I++) {  // Loop through all nodes where the line position and tension are to be computed
			if((s[I]<0.0) || (s[I]>L)) {
				Log->PrintfWarning("Warning from Catenary:: All line nodes must be located between the anchor and fairlead (inclusive) in routine Catenary(). Check ILine %d.",lineid);
				//cout << "        s[I]=" << s[I] << " and L=" << L << endl;
        return -1;
      }
      Ws=W*s[I]; //Initialize
      VFMinWLs=VFMinWL+Ws; // some commonly
      VFMinWLsOvrHF=VFMinWLs/HF; // used terms
      sOvrEA=s[I]/EA; // that depend
      SQRT1VFMinWLsOvrHF2=sqrt(1.0+VFMinWLsOvrHF*VFMinWLsOvrHF); // on s[I]

      if(s[I]<=LMinVFOvrW-HFOvrW/CB) { // .TRUE. if this node rests on the seabed and the tension is zero
        X[I]=s[I];
      	Z[I]=0.0;
				Te[I]=0.0;
			}
			else if(s[I]<=LMinVFOvrW) { // .TRUE. if this node rests on the seabed and the tension is nonzero
				X[I]=s[I]-(LMinVFOvrW-0.5*HFOvrW/CB)*HF/EA+sOvrEA*(HF+CB*VFMinWL+0.5*Ws*CB)+0.5*CB*VFMinWL*VFMinWL/WEA;
				Z[I]=0.0;
				Te[I]=HF+CB*VFMinWLs;
			}
			else { // LMinVFOvrW<s<=L ! This node must be above the seabed
				X[I]=log(VFMinWLsOvrHF+SQRT1VFMinWLsOvrHF2)*HFOvrW+sOvrEA*  HF+LMinVFOvrW-(LMinVFOvrW-0.5*HFOvrW/CB)*HF/EA;
				Z[I]=(-1.0+SQRT1VFMinWLsOvrHF2)*HFOvrW+sOvrEA*(VFMinWL+0.5*Ws)+0.5*   VFMinWL*VFMinWL/WEA;
				Te[I]=sqrt(HF*HF+VFMinWLs*VFMinWLs);
			}
		}
	}

	*HFout=HF;
	*VFout=VF;
	*HAout=HA;
	*VAout=VA;

	return 1;
}