(* ::Package:: *)

(* ::Title:: *)
(*Fundamental Mechanism Library*)


BeginPackage["fundMechanism`"]


(* ::Section:: *)
(*Kinematics*)


TraslMatrix::usage = "Matrice di traslazione"
Rotx::usage = "Matrice di rotazione"
Roty::usage= "Matrice di rotazione"
Rotz::usage= "Matrice di rotazione"

Begin["`Privat`"]
	TraslMatrix[x_,y_,z_]:=({
 {1, 0, 0, x},
 {0, 1, 0, y},
 {0, 0, 1, z},
 {0, 0, 0, 1}
});
	Rotx[\[Gamma]_]:=({
 {1, 0, 0, 0},
 {0, Cos[\[Gamma]], -Sin[\[Gamma]], 0},
 {0, Sin[\[Gamma]], Cos[\[Gamma]], 0},
 {0, 0, 0, 1}
}); 
	Roty[\[Gamma]_]:=({
 {Cos[\[Gamma]], 0, Sin[\[Gamma]], 0},
 {0, 1, 0, 0},
 {-Sin[\[Gamma]], 0, Cos[\[Gamma]], 0},
 {0, 0, 0, 1}
}); 
	Rotz[\[Gamma]_]:=({
 {Cos[\[Gamma]], -Sin[\[Gamma]], 0, 0},
 {Sin[\[Gamma]], Cos[\[Gamma]], 0, 0},
 {0, 0, 1, 0},
 {0, 0, 0, 1}
});
End[]


(* ::Section:: *)
(*Quadrilatero RRRR*)


(* ::Text:: *)
(*Si faccia riferimento a `AN.PVA_4L-4R.pdf`.*)


Quadrilatero::usage = "Quadrilatero[q,xA,yA,xD,yD,L1,L2,L3,modo], where A & D are pin joint in frame. It will return {\[Theta]2,\[Theta]3,ABCD}"

Begin["`Privat`"]
	Quadrilatero[q_,xA_,yA_,xD_,yD_,L1_,L2_,L3_,modo_]:= 
	(*q = movente, {\[Theta]2, \[Theta]3} = cedenti*)
		Module[
			{L5,xB,yB,xC,yC,\[Alpha],\[Theta]5,\[Theta]2,\[Theta]3, ABCD},
			(*1. Posizione di B*)
			xB = xA + L1*Cos[q];
			yB = yA + L1*Sin[q];
			(*2. Valore di L5 *)
			L5 =Sqrt[(xD-xB)^2 + (yD-yB)^2];
			(*3. T. del triangolo*)
			\[Alpha] = ArcCos[(L2^2+L5^2-L3^2)/(2*L2*L5)];
			(*4. Valore di \[Theta]5 -> \[Theta]2*)
			\[Theta]5 = ArcTan[xD-xB,yD-yB];
			\[Theta]2 = If[modo==-1,\[Theta]5-\[Alpha],\[Theta]5+\[Alpha]];
			(*5. Posizione di C*)
			xC=xB+L2*Cos[\[Theta]2];
			yC=yB+L2*Sin[\[Theta]2];
			(*6. Valore di \[Theta]3*)
			\[Theta]3=ArcTan[xC-xD,yC-yD];
			(*7. Ritorna il poligono e i cedenti*)
			ABCD = {{xA,yA},{xB,yB},{xC,yC},{xD,yD}};
			{\[Theta]2,\[Theta]3,ABCD}
		]
End[]


(* ::Section:: *)
(*Closing the package*)


EndPackage[]
