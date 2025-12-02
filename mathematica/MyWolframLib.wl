(* ::Package:: *)

(* ::Title:: *)
(*Wolfram Library*)


BeginPackage["MyWolframLib`"]


(* ::Section:: *)
(*Graphics functions*)


PinJointFrame::usage = "PinJointFrame[x,y,r] draws a pin joint in frame of center (x,y) and radius r"
PlotFrame::usage = "PlotFrame[T,l,o] plots the reference system's axes of absolute coordinates determinated by the traslation matrix T. The axes length = `l` = 5 as dafault. Axes's opacity = `o` = 0.7 as dafault."

(*$$ PinJointFrame $$*)
Begin["`Privat`"]
	PinJointFrame[x_,y_,r_]:=Graphics[
		{
			{Thick, Black, Circle[{x,y},r] (*inserire spessore eventualmente*)},
			{Black, Disk[{x,y},r,{0, Pi/2}]},
			{Black, Disk[{x,y},r,{Pi,3*Pi/2}]}
		}
	]
End[]

(*$$ PlotFrame $$*)
Begin["`Privat`"]
	PlotFrame[T_, l_:5, o_:0.7]:= Module[
		{O, X, Y, Z},
		O = T[[1;;3, 4]];
		X = O + l*T[[1;;3, 1]];
		Y = O + l*T[[1;;3, 2]];
		Z = O + l*T[[1;;3, 3]];
		Graphics3D[{ 
			Opacity[o],Thick, Red, Arrow[{O, X}],
			Opacity[o],Thick, Green, Arrow[{O, Y}],
			Opacity[o],Thick, Blue, Arrow[{O, Z}]
		}]
	]
End[]


(* ::Section:: *)
(*Closing the package*)


EndPackage[]
