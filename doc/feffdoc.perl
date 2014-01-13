## -*- mode: cperl -*-
## feffdoc.perl - handle feffdoc.sty in latex2html
## $Id: feffdoc.perl,v 1.1.1.1 2006/01/12 06:37:42 hebhop Exp $

package main;

$TITLE="FEFF8 Documentation";
$HTML_VERSION = 3.2;
$FEFF_VERSION = 8;
$FEFF_VNUM = 8.00;

## this handles the Card environment for the feff8 document.  the
## various substitutions handle the formatting of the arguments to the
## Card environment and the cross-referencing
sub do_env_Card {
  print "\n *** doing Card *** " if ($VERBOSITY > 1);
  local($_) = @_;
  local($anchor) = '';
  $_ = &translate_environments($_);
  $_ = &translate_commands($_);
  local($env_id) = " CLASS=\"CARD\"" if ($USING_STYLES);
  s!<\#\d+\#>([^<]*)<\#\d+\#>!<BIG><STRONG>$1</STRONG></BIG>&nbsp;&nbsp;&nbsp;\n!;
  s!<\#\d+\#>([^<]*)<\#\d+\#>!<KBD>$1</KBD>&nbsp;&nbsp;&nbsp;\n!;
  s!<\#\d+\#>([^<]*)<\#\d+\#>!<I>($1)</I><BR>\n!;
  (s!<\#\d+\#>([^<]*)<\#\d+\#>!!) && # set hypertext anchors
    ($anchor = &anchor_label("card:".$1,$CURRENT_FILE,''));
  $_ = $anchor . $_;
  $_;
}

## reference list environment -- this is a modified description list,
## so just use the html <DL>...</DL>
sub do_env_Reflist {
  print "\n *** doing Reflist *** " if ($VERBOSITY > 1);
  &do_env_description(@_);
}

## for some reason, the l2h slurps the abstract in with the \MakeTitle
## command, so I have to append $_ to the end of the title lines.
sub do_cmd_MakeTitle {
  local($_) = @_;
  my ($title_lines) =
    "<H1>FEFF$FEFF_VERSION</H1>\n" .
      "<BIG><STRONG>The FEFF Project<BR>\n" .
	"Department of Physics<BR>\n" .
	  "University of Washington<BR>\n" .
	    "Version $FEFF_VNUM</STRONG></BIG>\n" ;
  s!\\newchapter<\#\d+\#>([^<]*)<\#\d+\#>!!g;
  $title_lines . $_;
}

## set some commands to no-op in latex2html
sub do_cmd_newchapter {};
## sub do_cmd_tightlist {}; # no-op-ing this causes problems in Table 2.1
sub do_cmd_stretch {};

## I do not know how to place the module name on the description entry
## line, so I do this visually unappealing hack instead.
sub do_cmd_dotfill {
  local($_) = @_;
  "----------------" . $_;
};

## fonts
sub do_cmd_file {
  local($_) = @_;
  s!<\#\d+\#>([^<]*)<\#\d+\#>!<KBD>'$1'</KBD>!;
  $_;
}
sub do_cmd_command {
  local($_) = @_;
  s!<\#\d+\#>([^<]*)<\#\d+\#>!<KBD>$1</KBD>!;
  $_;
}
sub do_cmd_program {
  local($_) = @_;
  '<SMALL>' . &do_cmd_textsc($_) . '</SMALL>';
}
sub do_cmd_mathtt {
  local($_) = @_;
  &do_cmd_texttt($_);
}
sub do_cmd_textrm {
  local($_) = @_;
  $_;
}
sub do_cmd_mathrm {
  local($_) = @_;
  $_;
}
sub do_cmd_mathit {
  local($_) = @_;
  &do_cmd_textit($_);
}
sub do_cmd_mathcal {
  local($_) = @_;
  &do_cmd_textit($_);
}

## program names
sub do_cmd_feffcur {
  "<SMALL>FEFF</SMALL>" . $FEFF_VERSION;
}
sub do_cmd_feff {
  "<SMALL>FEFF</SMALL>";
}
sub do_cmd_feffit {
  "<SMALL>FEFFIT</SMALL>";
}
sub do_cmd_atoms {
  "<SMALL>ATOMS</SMALL>";
}


## math stuff
sub do_cmd_ell {
  local($_) = @_;
  "l";
};
sub do_cmd_star {
  local($_) = @_;
  "*";
};
sub do_cmd_cdot {
  local($_) = @_;
  "*";
};
sub do_cmd_Im {
  local($_) = @_;
  "Im" . $_;
};
sub do_cmd_Re {
  local($_) = @_;
  "Im" . $_;
};


1;
