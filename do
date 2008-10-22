#!/usr/bin/perl
# $Id$
#
# Usage: makemake {<program name> {<F90 compiler or fc or f77 or cc or c>}}
#
# Generate a Makefile from the sources in the current directory.  The source
# files may be in either C, FORTRAN 77, Fortran 90 or some combination of
# these languages.  If the F90 compiler specified is cray or parasoft, then
# the Makefile generated will conform to the conventions of these compilers.
# To run makemake, it will be necessary to modify the first line of this script
# to point to the actual location of Perl on your system.
#
# Written by Michael Wester <wester@math.unm.edu> February 16, 1995
# Cotopaxi (Consulting), Albuquerque, New Mexico
#
# Modified by Kacper Kowalik <kowalik@astri.umk.pl> October, 24, 2007
#
use Switch;
use File::Path;
use File::Copy;
my $argc = $#ARGV+1;
if($argc < 1){
   print "USAGE: makemake <problem>\n";
   die;
}
open MAKEIN, "< compilers/Makefile.in" or die $!;
@makein = <MAKEIN>;
close MAKEIN;
rmtree(['obj']);
mkpath(['obj']);
$probdir = "problems/". $ARGV[0] . "/";
if(!-r $probdir) { die "Can't open $probdir: $!";}

@prob = ( "../" . $probdir . "init_problem.F90", 
          "../" . $probdir . "problem.par" );
@base = <src/base/*.F90>;
for (@base) {
   s/src/\.\.\/src/;
}
@addons = ();
copy("compilers/newcompiler","obj/newcompiler");
system("chmod a+x obj/newcompiler");
$defs = $probdir."piernik.def";
copy($defs,"obj/piernik.def");
open (defs) or die "Can't open the file piernik.def!";
@fdefs = <defs>;
@d = grep (/define/,@fdefs);
if( grep { /GRAV/ }  @d) {
   push(@addons, "../src/gravity/gravity.F90");
   push(@addons, "../src/gravity/hydrostatic.F90");
   }
if( grep { /SELF_GRAV/ }  @d) {push(@addons, "../src/gravity/poisson_solver.F90");}
if( grep { /RESIST/} @d) {push(@addons, "../src/resist/resistivity.F90");}
if( grep { /SHEAR/}  @d) {push(@addons, "../src/shear/shear.F90");}
if( grep { /COSM_RAYS/} @d) {
   push(@addons, "../src/cosm_rays/cr_diffusion.F90");
}
if( grep { /SN_SRC/} @d) {
   push(@addons, "../src/supernovae/sn_sources.F90");
}
if( grep { /SNE_DISTR/} @d) {
   push(@addons, "../src/supernovae/sn_distr.F90");
}
if( grep { /ORIG/ } @d) {
   $solv="orig";
} else {
   $solv="ssp";
}
if( grep { /UNSPLIT_1D/ } @d) {
   push(@addons, "../src/scheme/unsplit/1D/mhdstep.F90");
   $sch="unsplit";
} elsif( grep { /UNSPLIT_3D/ } @d) {
   push(@addons, "../src/scheme/unsplit/3D/mhdstep.F90");
   $sch="unsplit";
} else {
   push(@addons, "../src/scheme/split/mhdstep.F90");
   $sch = 'split';
}
push(@addons, "../src/scheme/".$sch."/advects.F90");
push(@addons, "../src/scheme/".$sch."/fluids.F90");
push(@addons, "../src/scheme/".$sch."/".$solv."/tv.F90");
@files = ( @base, @prob, @addons );
@symln = ();
foreach $file (@files) {
   $pos = rindex($file,"\/")+1;
   $len = length($file);
   push(@symln,"obj/".substr($file,$pos,$len-$pos) );
}
for $i (0 .. $#symln){
  symlink($files[$i],$symln[$i]);
}
chdir 'obj';
open(MAKEFILE, "> Makefile-prep");

# Source listing
#
#print MAKEFILE @makein;
#foreach $line (@makein){
#   print "$line";
#}
print MAKEFILE "SRCS =\t";
@srcs = <*.F90>;
&PrintWords(8, 0, @srcs);
print MAKEFILE "\n\n";
#
# Object listing
#
print MAKEFILE "OBJS =\t";
@objs = @srcs;
foreach (@objs) { s/\.[^.]+$/.o/ };
&PrintWords(8, 0, @objs);
print MAKEFILE "\n\n";
#
# Define common macros
#
print MAKEFILE "LIBS = -L\$(HDF_LIB) -L\${MHDF_LIB} -lmfhdf -ldf -ljpeg -lz\n\n";
#
# make
#
print MAKEFILE "all: date \$(PROG) \n\n";
print MAKEFILE "\$(PROG): \$(OBJS)\n";
print MAKEFILE "\t\$(", &LanguageCompiler($ARGV[1], @srcs);
print MAKEFILE ") \$(LDFLAGS) -o \$@ \$(OBJS) \$(LIBS)\n\n";
#
# make date
#
print MAKEFILE "date: \n";
print MAKEFILE "\trm -f env.*\n";
print MAKEFILE "\tcp ../src/base/env.F90 .\n";
print MAKEFILE "\thead -n 1 *.F90 | grep Id > env.dat\n";
print MAKEFILE "\tcat piernik.def >> env.dat\n";
print MAKEFILE "\twc -l env.dat | awk '{print \"\tinteger, parameter :: nenv = \"\$\$1\"+1\"}' - >> env.F90\n";
print MAKEFILE "\techo \"   character*128, dimension(nenv) :: env = (/ &\" >> env.F90\n";
print MAKEFILE "\tawk '{printf(\"\\t\\t \\\" %s \\\" ,&\\n\",\$\$0)}' env.dat >> env.F90\n";
print MAKEFILE "\techo '      \" \" /)' >> env.F90\n";
print MAKEFILE "\techo 'end module comp_log' >> env.F90\n\n";

#
# make clean
#
print MAKEFILE "clean:\n";
print MAKEFILE "\trm -f \$(PROG) \$(OBJS) *.mod\n\n";
#
# make clean-run
#
print MAKEFILE "clean-run:\n";
print MAKEFILE "\trm -f *.bck *~ *.hdf *.res *.log *.tsl *.out *.tmp core*\n\n";
#
# make clean-all
#
print MAKEFILE "clean-all:\n";
print MAKEFILE "\trm -f \$(PROG) \$(OBJS) *.mod *.bck *~ *.hdf *.res *.log *.tsl *.out *.tmp core* *.f *.dbg \n\n";
#
# Make .F90 a valid suffix
#
print MAKEFILE ".SUFFIXES: \$(SUFFIXES) .F90\n\n";
#
# .F90 -> .o
#
print MAKEFILE ".F90.o:\n";
print MAKEFILE "\t\$(F90) \$(F90FLAGS) -c \$<\n\n";
#
# override the built-in rule for .mod (Modula-2 source code files)
#
print MAKEFILE "%.o : %.mod\n\n";
#
# Dependency listings
#
&MakeDependsf90($ARGV[1]);
&MakeDepends("*.f *.F *.F90", '^\s*use\s+["\']([^"\']+)["\']');
&MakeDepends("*.c",     '^\s*#\s*include\s+["\']([^"\']+)["\']');

system("./post","mhd",$ARGV[0]);
system("./newcompiler Makefile");

#
# &PrintWords(current output column, extra tab?, word list); --- print words
#    nicely
#
sub PrintWords {
   local($columns) = 78 - shift(@_);
   local($extratab) = shift(@_);
   local($wordlength);
   #
   print MAKEFILE @_[0];
   $columns -= length(shift(@_));
   foreach $word (@_) {
      $wordlength = length($word);
      if ($wordlength + 1 < $columns) {
         print MAKEFILE " $word";
         $columns -= $wordlength + 1;
         }
      else {
         #
         # Continue onto a new line
         #
         if ($extratab) {
            print MAKEFILE " \\\n\t\t$word";
            $columns = 62 - $wordlength;
            }
         else {
            print MAKEFILE " \\\n\t$word";
            $columns = 70 - $wordlength;
            }
         }
      }
   }

#
# &LanguageCompiler(compiler, sources); --- determine the correct language
#    compiler
#
sub LanguageCompiler {
   local($compiler) = &toLower(shift(@_));
   local(@srcs) = @_;
   #
   if (length($compiler) > 0) {
      CASE: {
         grep(/^$compiler$/, ("fc", "f77")) &&
            do { $compiler = "FC"; last CASE; };
         grep(/^$compiler$/, ("cc", "c"))   &&
            do { $compiler = "CC"; last CASE; };
         $compiler = "F90";
         }
      }
   else {
      CASE: {
         grep(/\.F90$/, @srcs)   && do { $compiler = "F90"; last CASE; };
         grep(/\.(f|F)$/, @srcs) && do { $compiler = "FC";  last CASE; };
         grep(/\.c$/, @srcs)     && do { $compiler = "CC";  last CASE; };
         $compiler = "???";
         }
      }
   $compiler;
   }

#
# &toLower(string); --- convert string into lower case
#
sub toLower {
   local($string) = @_[0];
   $string =~ tr/A-Z/a-z/;
   $string;
   }

#
# &uniq(sorted word list); --- remove adjacent duplicate words
#
sub uniq {
   local(@words);
   foreach $word (@_) {
      if ($word ne $words[$#words]) {
         push(@words, $word);
         }
      }
   @words;
   }

#
# &MakeDepends(language pattern, include file sed pattern); --- dependency
#    maker
#
sub MakeDepends {
   local(@incs);
   local($lang) = @_[0];
   local($pattern) = @_[1];
   #
   foreach $file (<${lang}>) {
      open(FILE, $file) || warn "Cannot open $file: $!\n";
      while (<FILE>) {
         /$pattern/i && push(@incs, $1);
         }
      if (defined @incs) {
         $file =~ s/\.[^.]+$/.o/;
         print MAKEFILE "$file: ";
         &PrintWords(length($file) + 2, 0, @incs);
         print MAKEFILE "\n";
         undef @incs;
         }
      }
   }

#
# &MakeDependsf90(f90 compiler); --- FORTRAN 90 dependency maker
#
sub MakeDependsf90 {
   local($compiler) = &toLower(@_[0]);
   local(@dependencies);
   local(%filename);
   local(@incs);
   local(@modules);
   local($objfile);
   #
   # Associate each module with the name of the file that contains it
   #
   foreach $file (<*.F90>) {
      open(FILE, $file) || warn "Cannot open $file: $!\n";
      while (<FILE>) {
         /^\s*module\s+([^\s!]+)/i &&
            ($filename{&toLower($1)} = $file) =~ s/\.F90$/.o/;
         }
      }
   #
   # Print the dependencies of each file that has one or more include's or
   # references one or more modules
   #
   foreach $file (<*.F90>) {
      open(FILE, $file);
      while (<FILE>) {
#         /^\s*include\s+["\']([^"\']+)["\']/i && push(@incs, $1);
         /^\s*use\s+([^\s,!]+)/i && push(@modules, &toLower($1));
         }
      if (defined @incs || defined @modules) {
         ($objfile = $file) =~ s/\.F90$/.o/;
         print MAKEFILE "$objfile: $file ";
         undef @dependencies;
         foreach $module (@modules) {
            push(@dependencies, $filename{$module});
            }
         @dependencies = &uniq(sort(@dependencies));
         &PrintWords(length($objfile) + 2, 0,
                     @dependencies, &uniq(sort(@incs)));
         print MAKEFILE "\n";
         undef @incs;
         undef @modules;
         #
         # Cray F90 compiler
         #
         if ($compiler eq "cray") {
            print MAKEFILE "\t\$(F90) \$(F90FLAGS) -c ";
            foreach $depend (@dependencies) {
               push(@modules, "-p", $depend);
               }
            push(@modules, $file);
            &PrintWords(30, 1, @modules);
            print MAKEFILE "\n";
            undef @modules;
            }
         #
         # ParaSoft F90 compiler
         #
         if ($compiler eq "parasoft") {
            print MAKEFILE "\t\$(F90) \$(F90FLAGS) -c ";
            foreach $depend (@dependencies) {
               $depend =~ s/\.o$/.F90/;
               push(@modules, "-module", $depend);
               }
            push(@modules, $file);
            &PrintWords(30, 1, @modules);
            print MAKEFILE "\n";
            undef @modules;
            }
         }
      }
   }
