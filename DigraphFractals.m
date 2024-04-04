(* :Name: DigraphFractals` *)
(* :Author: Mark McClure, November 1999 *)
(* :Copyright: Copyright Mark McClure, 1999 *)
(* :Package Version: 1.0 *)
(* :Mathematica Version: 4.0 *)
(* :Summary:
	Generates images of Digraph Self Similar sets. 
*)

BeginPackage["DigraphFractals`"]

Needs["Utilities`FilterOptions`"]

DigraphFractals::usage = "DigraphFractals is a package defining
several functions used to generate images of digraph Self 
Similar sets."

ShowDigraphFractals::usage = "ShowDigraphFractals[digraph_,depth_]
generates the approximation to the digraph Self Similar sets 
defined by digraph to order depth using a deterministic algorithm.
The digraph is represented by a matrix of lists of affine functions
written in the form {{{a1,a2},{b1,b2}},{x0,y0}}."

ShowDigraphFractalsStochastic::usage = 
"ShowDigraphFractalsStochastic[digraph_,numPoints_]
generates numPoints points approximating the digraph Self Similar
sets defined by digraph using a stochastic algorithm.
The digraph is represented by a matrix of lists of affine functions
written in the form {{{a1,a2},{b1,b2}},{x0,y0}}."

ComputePMatrix::usage = "ComputePMatrix[digraph_] estimates a 
probability matrix for digraph which generates a uniform distribution
of points over each digraph Fractal."

ComputeDimension::usage = "ComputeDimension[digraph_] computes a 
numerical estimate of the dimension of the corresponding digraph
fractals."

Initiators::usage = "Initiators is an option for ShowDigraphFractals
indicating an initial list of Graphics primitives for the digraph
to operate on.  Initiators should be a list of lists of Graphics
primitives with length equal to the dimension of the digraph."

PMatrix::usage = "PMatrix is an option for 
ShowDigraphFractalsStochastic indicating a matrix of probabilities
for the functions in the digraph. in
ShowDigraphFractalsStochastic[digraph_, numPoints_], PMatrix->pMatrix],
pMatrix should be a matrix of lists of positive numbers with the same 
shape as digraph.  The sum of the sums of the lists in any column of 
PMatrix should be one."

Begin["`Private`"]

(* Set the Options *)

(* For ShowDigraphFractals *)
ShowDigraphFractalsOptions = Join[{AspectRatio -> Automatic,
		Initiators -> Automatic, Prolog -> AbsolutePointSize[.05]},
		Options[Graphics]];
keywords = Complement[Union[First /@ ShowDigraphFractalsOptions], 
	{DisplayFunction, DefaultFont}];
vals = keywords /.  ShowDigraphFractalsOptions;
special = {DisplayFunction :> $DisplayFunction, DefaultFont :> $DefaultFont};
Options[ShowDigraphFractals] = 
	Union[Apply[Rule,Transpose[{keywords,vals}],{1}], special];

(* For ShowDigraphFractalsStochastic *)
ShowDigraphFractalsStochasticOptions = Join[{AspectRatio -> Automatic,
	Prolog -> AbsolutePointSize[.05], PMatrix -> Automatic}, Options[Graphics]];
keywords = Complement[Union[First /@ ShowDigraphFractalsStochasticOptions], 
	{DisplayFunction, DefaultFont}];
vals = keywords /.  ShowDigraphFractalsStochasticOptions;
special = {DisplayFunction :> $DisplayFunction, DefaultFont :> $DefaultFont};
Options[ShowDigraphFractalsStochastic] = 
	Union[Apply[Rule,Transpose[{keywords,vals}],{1}], special];


(* The Functions *)

(* Error Messages *)
DigraphFractals::optx = "Unknown Option `1`."
DigraphFractals::badDigraph = "The first argument in `1` must be a
matrix of lists of affine lists."
DigraphFractals::badInt = "The second argument in `1` must be a
non-negative integer."
DigraphFractals::badOpt = "Options expected as optional arguments in `1`."


ShowDigraphFractals[digraph_, depth_, opts___] :=
	
	Module[{digraphFuncs, toFuncs, placeHolder, start, 
		graphicList, initiators, temp,
		valid = First /@ Options[ShowDigraphFractals]},

	Scan[If[!MemberQ[valid, First[#]],
		Message[DigraphFractals::optx, ToString[First[#]]]]&,
		Flatten[{opts}]
	];

	initiators = Initiators /. {opts} /. Options[ShowDigraphFractals];
	If[initiators === Automatic,
		initiators = Table[{Point[{0,0}]}, {Length[digraph]}]
	];
	

	toFuncs[{A_,b_}] := Module[{fOut},
		fOut[{x_,y_}] := N[A.{x,y} + b] /; NumericQ[x] && NumericQ[y];
		fOut[Point[{x_,y_}]] := Point[fOut[{x,y}]];
		fOut[x_List] := fOut /@ x;
		fOut[Line[x__]] := Line[fOut /@ x];
		fOut[Polygon[x__]] := Polygon[fOut /@ x];
		fOut[Disk[p_,r_]] := Disk[fOut[p],r];
		fOut[Circle[p_,r_]] := Circle[fOut[p],r];
		
		(* Why Doesn't This Work !??! *)
		fOut[Arrow[arrowStart_, arrowFinish_, arrowOpts___]] := 
			Arrow[fOut[arrowStart], fOut[arrowFinish], arrowOpts];
		
		fOut[h_[x__]] := h[x];
		fOut
	];
		
	start = placeHolder @@ # & /@ initiators;
	digraphFuncs = Map[toFuncs,digraph, {3}];
	digraphFuncs = Map[placeHolder @@ #  &, digraphFuncs, {2}] ;
	graphicList = Nest[Inner[Outer[#1[#2]&, #1, #2]&, digraphFuncs, #,List]&,
			start, depth] //. placeHolder[x___] -> {x};
	Show[Graphics[#], FilterOptions[Graphics,opts],
		FilterOptions[Graphics, 
			Sequence @@ Options[ShowDigraphFractals]]]& /@ graphicList
] /;
	(And @@ Flatten[Map[MatchQ[#,
		{{{a_?NumericQ, b_?NumericQ}, {c_?NumericQ, d_?NumericQ}},
		{e_?NumericQ, f_?NumericQ}}] &, digraph, {3}]] ||
		Message[DigraphFractals::badDigraph, ShowDigraphFractals] )  &&
	( (IntegerQ[depth] && depth >= 0 ) ||
		Message[DigraphFractals::badInt, ShowDigraphFractals, depth] ) &&
	( And @@ Map[OptionQ, {opts}] ||
	  	Message[DigraphFractals::badOpt,ShowDigraphFractals] )


ShowDigraphFractalsStochastic[digraph_, numPoints_, opts___] :=
	
  Module[{matrices, eigenvalues, s, sMatrix,
		spectralRadius, spectralRadius0,
		approximateDim, dimMatrix, perronNumbers,
		pMatrix1, pMatrix2, pMatrixNormalizer,
		toFuncs, digraphFuncs, v, v1, v2,
		currentPoint, points, choose,
		valid = First /@ Options[ShowDigraphFractalsStochastic]},

	Scan[If[!MemberQ[valid, First[#]],
		Message[DigraphFractals::optx, ToString[First[#]]]]&,
		Flatten[{opts}]
	];

	(* Computing the Probabilities *)
	
	pMatrix = PMatrix /. {opts} /. Options[ShowDigraphFractalsStochastic];
	
	pMatrixNormalizer[subPList_List] := If[
		Length[subPList] > 0,
			FoldList[Plus, 0, subPList]/Plus @@ subPList,
			subPList];
	If[pMatrix === Automatic,
	    {matrices = Map[First, Transpose[digraph], {3}] // N;
		eigenMatrix = matrices /. 
			{{a_Real, b_Real}, {c_Real, d_Real}} -> 
			Max[Abs /@ Eigenvalues[{{a, b}, {c, d}}]]^s;
		sMatrix = Map[Plus @@ # &, eigenMatrix, {2}];
		eigenvalues = Eigenvalues[sMatrix];
		spectralRadius0 = Max[Chop[Abs@N[eigenvalues /. s -> 0]]];
		spectralRadius = Select[eigenvalues, 
			((# /. s -> 0) == spectralRadius0) &][[1]];
		approximateDim = Chop[s /.
			FindRoot[spectralRadius == 1, {s, 1}]];
		dimMatrix = sMatrix /. s -> approximateDim;
		perronNumbers = Eigensystem[dimMatrix][[2]][[1]];
		pMatrix1 = Inner[Times, dimMatrix,
			perronNumbers, List]/perronNumbers;
		pMatrix1 = FoldList[Plus, 0, #] & /@ pMatrix1;
		pMatrix2 = eigenMatrix /. s -> approximateDim;
		pMatrix2 = Map[pMatrixNormalizer, pMatrix2, {2}];},
		
		{pMatrix = N[pMatrix];
		pMatrix2 = Transpose[Map[pMatrixNormalizer, pMatrix, {2}]];
		pMatrix1 = Map[Plus @@ # & , Transpose[pMatrix], {2}];
		pMatrix1 = FoldList[Plus, 0, #] & /@ pMatrix1;}
	];

	toFuncs[{A_, b_}] := A.# + b &;
	digraphFuncs = Map[toFuncs, N[digraph], {3}];
	v = 1; v1 = pMatrix1[[v]]; v2 = pMatrix2[[v]];
	currentPoint = {0., 0.};
	points = Table[{}, {Length[digraph]}];
	choose[v1_, v2_] := Module[
		{i, j, chooser},
		chooser = Random[];
		i = Length[Select[v1, (# < chooser) &]];
		chooser = Random[];
		j = Length[Select[v2[[i]], (# < chooser) &]];
		{i, j}];

	Do[{
		{i, j} = choose[v1, v2];
		currentPoint = 
			digraphFuncs[[i, v, j]][currentPoint];
		points[[i]] = {points[[i]], currentPoint};
		v = i;  v1 = pMatrix1[[v]];  v2 = pMatrix2[[v]];
		}, {numPoints}];
	
	points = Partition[Flatten[#], 2] & /@ points;
	Show[Graphics[Point /@ #], FilterOptions[Graphics,opts],
		FilterOptions[Graphics, 
			Sequence @@ Options[ShowDigraphFractalsStochastic]]] & /@ points 
] /;
	(And @@ Flatten[Map[MatchQ[#,
		{{{a_?NumericQ, b_?NumericQ}, {c_?NumericQ, d_?NumericQ}},
		{e_?NumericQ, f_?NumericQ}}] &, digraph, {3}]] ||
		Message[DigraphFractals::badDigraph, 
			ShowDigraphFractalsStochastic] )  &&
	( (IntegerQ[numPoints] && numPoints >= 0 ) ||
		Message[DigraphFractals::badInt, 
			ShowDigraphFractalsStochastic, numPoints] ) &&
	( And @@ Map[OptionQ, {opts}] ||
	  	Message[DigraphFractals::badOpt,ShowDigraphFractalsStochastic] )


ComputePMatrix[digraph_] := 

	Module[{
		matrices, eigenvalues, s, sMatrix,
		spectralRadius, spectralRadius0,
		approximateDim, dimMatrix, perronNumbers,
		pMatrix1, pMatrix2, pMatrixNormalizer},
		
		pMatrixNormalizer[subPList_List] := If[
			Length[subPList] > 0,
				subPList/Plus @@ subPList,subPList];
		matrices = Map[First, Transpose[digraph], {3}] // N;
		eigenMatrix = matrices /. 
			{{a_Real, b_Real}, {c_Real, d_Real}} -> 
			Max[Abs /@ Eigenvalues[{{a, b}, {c, d}}]]^s;
		sMatrix = Map[Plus @@ # &, eigenMatrix, {2}];
		eigenvalues = Eigenvalues[sMatrix];
		spectralRadius0 = Max[Chop[N[eigenvalues /. s -> 0]]];
		spectralRadius = Select[eigenvalues, 
			((# /. s -> 0) == spectralRadius0) &][[1]];
		approximateDim = Chop[s /.
			FindRoot[spectralRadius == 1, {s, 1}]];
		dimMatrix = sMatrix /. s -> approximateDim;
		perronNumbers = Eigensystem[dimMatrix][[2]][[1]];
		pMatrix1 = Inner[Times, dimMatrix,
			perronNumbers, List]/perronNumbers;
		pMatrix2 = eigenMatrix /. s -> approximateDim;
		pMatrix2 = Map[pMatrixNormalizer, pMatrix2, {2}];
		Transpose[pMatrix1 pMatrix2]
] /;
	(And @@ Flatten[Map[MatchQ[#,
		{{{a_?NumericQ, b_?NumericQ}, {c_?NumericQ, d_?NumericQ}},
		{e_?NumericQ, f_?NumericQ}}] &, digraph, {3}]] ||
		Message[DigraphFractals::badDigraph, 
			ComputePMatrix] )


ComputeDimension[digraph_] := Module[
	{matrices, eigenvalues, s, sMatrix,
		spectralRadius, spectralRadius0,
		approximateDim},

	matrices = Map[First, Transpose[digraph], {3}] // N;
	eigenMatrix = matrices /. 
		{{a_Real, b_Real}, {c_Real, d_Real}} -> 
		Max[Abs /@ Eigenvalues[{{a, b}, {c, d}}]]^s;
	sMatrix = Map[Plus @@ # &, eigenMatrix, {2}];
	eigenvalues = Eigenvalues[sMatrix];
	spectralRadius0 = Max[Chop[Abs@N[eigenvalues /. s -> 0]]];
	spectralRadius = Select[eigenvalues, 
		((# /. s -> 0) == spectralRadius0) &][[1]];
	approximateDim = Chop[s /.
		FindRoot[spectralRadius == 1, {s, 1}]]
] /;
	(And @@ Flatten[Map[MatchQ[#,
		{{{a_?NumericQ, b_?NumericQ}, {c_?NumericQ, d_?NumericQ}},
		{e_?NumericQ, f_?NumericQ}}] &, digraph, {3}]] ||
		Message[DigraphFractals::badDigraph, 
			ComputeDimension] )


End[]  (* End Private Context *)

Protect[DigraphFractals, ShowDigraphFractals, 
	ShowDigraphFractalsStochastic, ComputePMatrix, 
	Initiators, PMatrix]

EndPackage[]

(*
DigraphFractals`
*)
(*
DigraphFractals is a package defining several functions used to generate images of digraph Self  Similar sets.
\
*)
(*
ShowDigraphFractals[digraph_,depth_] generates the approximation to the digraph Self Similar sets  defined by digraph to order depth using a deterministic algorithm. The digraph is represented by a matrix of lists of affine functions written in the form {{{a1,a2},{b1,b2}},{x0,y0}}.
\
*)
(*
ShowDigraphFractalsStochastic[digraph_,numPoints_] generates numPoints points approximating the digraph Self Similar sets defined by digraph using a stochastic algorithm. The digraph is represented by a matrix of lists of affine functions written in the form {{{a1,a2},{b1,b2}},{x0,y0}}.
\
*)
(*
ComputePMatrix[digraph_] estimates a  probability matrix for digraph which generates a uniform distribution of points over each digraph Fractal.
\
*)
(*
ComputeDimension[digraph_] computes a  numerical estimate of the dimension of the corresponding digraph fractals.
\
*)
(*
Initiators is an option for ShowDigraphFractals indicating an initial list of Graphics primitives for the digraph to operate on.  Initiators should be a list of lists of Graphics primitives with length equal to the dimension of the digraph.
\
*)
(*
PMatrix is an option for  ShowDigraphFractalsStochastic indicating a matrix of probabilities for the functions in the digraph. in ShowDigraphFractalsStochastic[digraph_, numPoints_], PMatrix->pMatrix], pMatrix should be a matrix of lists of positive numbers with the same  shape as digraph.  The sum of the sums of the lists in any column of  PMatrix should be one.
\
*)
(*
DigraphFractals`Private`
*)
(*
Unknown Option `1`.
*)
(*
The first argument in `1` must be a matrix of lists of affine lists.
*)
(*
The second argument in `1` must be a non-negative integer.
*)
(*
Options expected as optional arguments in `1`.
*)
(*
\!\(\*
  RowBox[{\(General::"spell1"\), \(\(:\)\(\ \)\), "\<\"Possible spelling error: new symbol name \\\"\\!\\(initiators\\)\\\" is similar to existing symbol \\\"\\!\\(Initiators\\)\\\". \\!\\(\\*ButtonBox[\\\"More\[Ellipsis]\\\", ButtonStyle->\\\"RefGuideLinkText\\\", ButtonFrame->None, ButtonData:>\\\"General::spell1\\\"]\\)\"\>"}]\)
*)
(*
\!\(\*
  RowBox[{\(General::"spell1"\), \(\(:\)\(\ \)\), "\<\"Possible spelling error: new symbol name \\\"\\!\\(fOut\\)\\\" is similar to existing symbol \\\"\\!\\(Out\\)\\\". \\!\\(\\*ButtonBox[\\\"More\[Ellipsis]\\\", ButtonStyle->\\\"RefGuideLinkText\\\", ButtonFrame->None, ButtonData:>\\\"General::spell1\\\"]\\)\"\>"}]\)
*)
(*
\!\(\*
  RowBox[{\(General::"spell1"\), \(\(:\)\(\ \)\), "\<\"Possible spelling error: new symbol name \\\"\\!\\(eigenvalues\\)\\\" is similar to existing symbol \\\"\\!\\(Eigenvalues\\)\\\". \\!\\(\\*ButtonBox[\\\"More\[Ellipsis]\\\", ButtonStyle->\\\"RefGuideLinkText\\\", ButtonFrame->None, ButtonData:>\\\"General::spell1\\\"]\\)\"\>"}]\)
*)
(*
\!\(\*
  RowBox[{\(General::"spell"\), \(\(:\)\(\ \)\), "\<\"Possible spelling error: new symbol name \\\"\\!\\(pMatrix\\)\\\" is similar to existing symbols \\!\\({PMatrix, sMatrix}\\). \\!\\(\\*ButtonBox[\\\"More\[Ellipsis]\\\", ButtonStyle->\\\"RefGuideLinkText\\\", ButtonFrame->None, ButtonData:>\\\"General::spell\\\"]\\)\"\>"}]\)
*)
(*
\!\(\*
  RowBox[{\(General::"spell1"\), \(\(:\)\(\ \)\), "\<\"Possible spelling error: new symbol name \\\"\\!\\(chooser\\)\\\" is similar to existing symbol \\\"\\!\\(choose\\)\\\". \\!\\(\\*ButtonBox[\\\"More\[Ellipsis]\\\", ButtonStyle->\\\"RefGuideLinkText\\\", ButtonFrame->None, ButtonData:>\\\"General::spell1\\\"]\\)\"\>"}]\)
*)
(*
DigraphFractals`Private`
*)
(*
{DigraphFractals, ShowDigraphFractals, ShowDigraphFractalsStochastic, \
ComputePMatrix, Initiators, PMatrix}
*)
