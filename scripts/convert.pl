#!/usr/bin/perl -w
#
# script to convert files from old naming schemes to new ones
#


# names that are probably OK, but could do with removal
# of the _conv prefix
@removeconvprefix = (
"InitGridDef",
"AllocGridQuant",
"InitGridQuant",
"InitGridQuantSub",
"InitGridQuantLHAPDF",
"PrintGridQuant",
"EvalGridQuant",
"MomGridQuant",
"WgtGridQuant",
"AllocGridConv",
"InitGridConv",
"ValidateGD",
"GetDerivedProbes",
"SetDerivedConv",
"SetDerivedConv_nodealloc",
);

@removepdfgenprefix = (
"AllocPDF",
"InitPDF",
"InitPDFSub",
"InitPDF_LHAPDF",
"AllocInitPDF",
"AllocInitPDFSub",
);

# names that need to made more sensible (in convolution.f90)
$rename{qr/conv\\?_DelGridQuant/} = "Delete";
$rename{qr/conv\\?_DelGridConv/}  = "Delete";
$rename{qr/conv\\?_ZeroGridConv/} = "SetToZero";
$rename{qr/conv\\?_MultGridConv/} = "Multiply";
$rename{qr/conv\\?_AddGridConv/}  = "AddWithCoeff";
$rename{qr/conv\\?_ConvGridConv/} = "SetToConvolution";
$rename{qr/conv\\?_CommGridConv/} = "SetToCommutator";
$rename{qr/conv\\?_Seteps/} = "SetConvolutionEps";



# names that need to made more sensible (in conv_objects)
$rename{qr/cobj\\?_InitSplitPolLO/}   = "InitSplitMatPolLO";
$rename{qr/cobj\\?_InitSplitPolNLO/}  = "InitSplitMatPolNLO";
$rename{qr/cobj\\?_InitSplitLO/}   = "InitSplitMatLO";
$rename{qr/cobj\\?_InitSplitNLO/}  = "InitSplitMatNLO";
$rename{qr/cobj\\?_InitSplitNNLO/} = "InitSplitMatNNLO";
$rename{qr/cobj\\?_InitMTMNNLO/} = "InitMTMNNLO";

$rename{qr/cobj\\\?_GetDerivedProbes/} = "GetDerivedSplitMatProbes";

$result = '';

if ($#ARGV >= 0) {
  open ($filehandle, "<$ARGV[0]") || die "Could not open $ARGV[0] for reading";
} else {
  $filehandle = STDIN;
}

while ($line=<$filehandle>) {
  $result .= $line;
}

foreach $name (@removeconvprefix) {
  $prefixname = qr/conv\\?_$name/;
  $result =~ s/\b$prefixname\b/$name/gi;
}
foreach $name (@removepdfgenprefix) {
  $prefixname = qr/pdfgen\\?_$name/;
  $result =~ s/\b$prefixname\b/$name/gi;
}

foreach $name (keys %rename) {
  $newname = $rename{$name};
  $result =~ s/\b$name\b/$newname/gi;
}

# some more delicate cases
$result =~ s/\bid_([a-z]+)\b/iflv_$1/gi;



close($filehandle);

if ($#ARGV >= 0) {
  $name = $ARGV[0];
  print STDOUT "Renaming $name to $name.bak\n";
  rename ($name, "$name.bak");
  open (OUT, ">$name") || die "Could not open $name";
  print STDOUT "Writing new version to $name\n";
  print OUT $result;
} else {
  print $result;
}
