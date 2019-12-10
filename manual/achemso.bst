%%
%% This is file `achemso.bst',
  %% generated with the docstrip utility.
  %%
  %% The original source files were:
  %%
  %% achemso.dtx  (with options: `huvudkod')
  %% 
  %% Copyright (c) 1995-98 by Mats Dahlgren .
  %% All rights reserved.  See the file `achemso.ins' for information
  %% on how you may (re-)distribute the `achemso' package files.
  %% You are not allowed to make any changes to this file without
  %% explicit permission from the author.
  %% 
  ENTRY
    { address
          author
	        booktitle
		    chapter
		        edition
			    editor
			        howpublished
				    institution
				        journal
					    key
					        month
						    note
						        number
							    organization
							        pages
								    publisher
								        school
									    series
									        title
										    type
										        volume
											    year
											        annotate
												  }
  {}
    { label }

    INTEGERS { output.state before.all mid.sentence after.sentence after.block }

    FUNCTION {init.state.consts}
{ #0 'before.all :=
    #1 'mid.sentence :=
    #2 'after.sentence :=
    #3 'after.block :=
}

STRINGS { s t }

FUNCTION {output.nonnull}
{ 's :=
    output.state mid.sentence =
          { " " * write$ }
        { output.state after.block =
	          { ", " * write$
		         newline$
			        "\newblock " write$
				   }
		     { output.state before.all =
		              'write$
				       { add.period$ " " * write$ }
			           if$
				        }
		           if$
			           mid.sentence 'output.state :=
				       }
	  if$
	      s
}

FUNCTION {output}
{ duplicate$ empty$
      'pop$
	    'output.nonnull
	      if$
}

FUNCTION {output.check}
{ 't :=
    duplicate$ empty$
          { pop$ "empty " t * " in " * cite$ * warning$ }
        'output.nonnull
	    if$
}

FUNCTION {output.bibitem}
{ newline$
    "\bibitem{" write$
        cite$ write$
	  "}" write$
	    newline$
	      ""
	        before.all 'output.state :=
}

FUNCTION {fin.entry}
{ add.period$
    write$
        newline$
}

FUNCTION {new.block}
{ output.state before.all =
      'skip$
	    { after.block 'output.state := }
        if$
}

FUNCTION {new.sentence}
{ skip$
}

FUNCTION {not}
{   { #0 }
      { #1 }
        if$
}

FUNCTION {and}
{   'skip$
      { pop$ #0 }
        if$
}

FUNCTION {or}
{   { pop$ #1 }
      'skip$
	  if$
}

FUNCTION {new.block.checka}
{ empty$
      'skip$
	    'new.block
	      if$
}

FUNCTION {new.block.checkb}
{ empty$
    swap$ empty$
        and
	    'skip$
	        'new.block
		  if$
}

FUNCTION {field.or.null}
{ duplicate$ empty$
      { pop$ "" }
          'skip$
	      if$
}

FUNCTION {emphasize}
{ duplicate$ empty$
      { pop$ "" }
          { "\textit{" swap$ * "}" * }
	    if$
}

FUNCTION {emphasizec}
{ duplicate$ empty$
      { pop$ "" }
          { "\textit{" swap$ * ",}" * }
	    if$
}

FUNCTION {emphasizesc}
{ duplicate$ empty$
      { pop$ "" }
          { "\textit{" swap$ * ";}" * }
	    if$
}

FUNCTION {slantc}
{ duplicate$ empty$
      { pop$ "" }
          { "\textsl{" swap$ * ",}" * }
	    if$
}

FUNCTION {emboldenc}
{ duplicate$ empty$
      { pop$ "" }
          { "\textbf{" swap$ * ",}" * }
	    if$
}

INTEGERS { nameptr namesleft numnames }

INTEGERS { etal }

FUNCTION {format.names}
{ 's :=
    #1 'nameptr :=
    s num.names$ 'numnames :=
        numnames #15 >
	  s numnames "{ll}" format.name$ "others" = numnames #1 > and
	    or 'etal :=
	      etal
	          { #1 #1 + 'namesleft := }
        { numnames 'namesleft := }
	  if$
	        { namesleft #0 > }
	      { s nameptr "{vv~}{ll,}{~jj,}{~f.}" format.name$ 't :=
		      nameptr #1 >
			   { namesleft #1 >
			            { ";\ \ " * t * }
				           { nameptr #2 >
					           { "" * }
						         'skip$
							            if$
								             t "others" =
									              etal or
										            { " \textit{et al.}\ " * }
							       { ";\ \ " * t * }
							                if$
									         }
					        if$
						     }
		         't
			         if$
				       nameptr #1 + 'nameptr :=
				             namesleft #1 - 'namesleft :=
					         }
	        while$
}

FUNCTION {format.authors}
{ author empty$
      { "" }
          { author format.names }
	    if$
}

FUNCTION {format.editors}
{ editor empty$
      { "" }
          { editor format.names
	          editor num.names$ #1 >
		       { ",\ \ Eds.; " * }
		     { ",\ \ Ed.; " * }
		           if$
			         }
	    if$
}

FUNCTION {n.dashify}
{ 't :=
    ""
          { t empty$ not }
        { t #1 #1 substring$ "-" =
	     { t #1 #2 substring$ "--" = not
	              { "--" *
			         t #2 global.max$ substring$ 't :=
				          }
		             {   { t #1 #1 substring$ "-" = }
			             { "-" *
				               t #2 global.max$ substring$ 't :=
						       }
				              while$
						       }
			          if$
				       }
	        { t #1 #1 substring$ *
		       t #2 global.max$ substring$ 't :=
			    }
		      if$
			    }
	  while$
}

FUNCTION {format.date}
{ year empty$
      { "" }
          'year
	      if$
}

FUNCTION {tie.or.space.connect}
{ duplicate$ text.length$ #3 <
      { "~" }
          { " " }
	    if$
	        swap$ * *
}

FUNCTION {either.or.check}
{ empty$
      'pop$
	    { "can't use both " swap$ * " fields in " * cite$ * warning$ }
        if$
}

FUNCTION {format.bvolume}
{ volume empty$
      { "" }
          { "volume" volume tie.or.space.connect
	          series empty$
		       'skip$
		          { " of " * series emphasize * }
		        if$
			        "volume and number" number either.or.check
				    }
	    if$
}

FUNCTION {format.number.series}
{ volume empty$
      { number empty$
	   { series field.or.null }
	           { series empty$
		            { "there's a number but no series in " cite$ * warning$ }
			           { series " " * number * "; " * }
				        if$
					     }
		         if$
			       }
          { "" }
	    if$
}

FUNCTION {format.edition}
{ edition empty$
      { "" }
          { edition " ed.;" * }
	    if$
}

INTEGERS { multiresult }

FUNCTION {multi.page.check}
{ 't :=
    #0 'multiresult :=
      { multiresult not
	      t empty$ not
		      and
		          }
          { t #1 #1 substring$
	          duplicate$ "-" =
		          swap$ duplicate$ "," =
			        swap$ "+" =
				      or or
				         { #1 'multiresult := }
		     { t #2 global.max$ substring$ 't := }
		           if$
			         }
	    while$
	        multiresult
}

FUNCTION {format.pages}
{ pages empty$
      { "" }
          { pages multi.page.check
	       { "pages" pages n.dashify tie.or.space.connect }
	          { "page" pages tie.or.space.connect }
		        if$
			      }
	    if$
}

FUNCTION {format.pages.a}
{ pages empty$
      { "" }
          { "page" pages tie.or.space.connect }
	    if$
}

FUNCTION {format.yr.vol.pages}
{ year field.or.null emboldenc
    " " swap$ * *
        volume field.or.null slantc
	  " " swap$ * *
	    pages empty$
	        'skip$
		    { duplicate$ empty$
		         { pop$ format.pages.a }
			    { " " * pages * }
			          if$
				        }
      if$
}

FUNCTION {format.chapter.pages}
{ chapter empty$
      'format.pages
	    { type empty$
	         { "Chapter" }
		    { type "l" change.case$ }
		          if$
			          chapter tie.or.space.connect
				        pages empty$
					   'skip$
					      { ", " * format.pages * }
			        if$
				      }
        if$
}

FUNCTION {empty.misc.check}
{ author empty$ title empty$ howpublished empty$
    month empty$ year empty$ note empty$
        and and and and and
	    { "all relevant fields are empty in " cite$ * warning$ }
        'skip$
	    if$
}

FUNCTION {format.thesis}
{ type empty$
      {"Thesis,"}
          {type " thesis," * }
	    if$
}

FUNCTION {format.article.crossref}
{ key empty$
      { journal empty$
	   { "need key or journal for " cite$ * " to crossref " * crossref *
	          warning$
		         ""
			    }
	      { "In " journal * }
	            if$
		          }
          { "In " key * }
	    if$
	        " \cite{" * crossref * "}" *
}

FUNCTION {format.crossref.editor}
{ editor #1 "{vv~}{ll}" format.name$
    editor num.names$ duplicate$
        #2 >
          { pop$ " et~al." * }
        { #2 <
	     'skip$
	          { editor #2 "{ff }{vv }{ll}{ jj}" format.name$ "others" =
		           { " et~al." * }
			          { " and " * editor #2 "{vv~}{ll}" format.name$ * }
				       if$
					    }
	           if$
		         }
	  if$
}

FUNCTION {format.book.crossref}
{ volume empty$
      { "empty volume in " cite$ * "'s crossref of " * crossref * warning$
	      "In "
		    }
          { "Volume" volume tie.or.space.connect
	          " of " *
		        }
	    if$
	        editor empty$
		  editor field.or.null author field.or.null =
		    or
		        { key empty$
			     { series empty$
			              { "need editor, key, or series for " cite$ * " to crossref " *
					         crossref * warning$
						            "" *
							           }
				             { "\textit{" * series * "}" * }
					          if$
						       }
			        { key * }
				      if$
					    }
	        { format.crossref.editor * }
		  if$
		      " \cite{" * crossref * "}" *
}

FUNCTION {format.incoll.inproc.crossref}
{ editor empty$
    editor field.or.null author field.or.null =
        or
	    { key empty$
	         { booktitle empty$
		          { "need editor, key, or booktitle for " cite$ * " to crossref " *
			             crossref * warning$
				                ""
						       }
			         { "In \textit{" booktitle * ",}" * }
				      if$
					   }
		    { "In " key * }
		          if$
			        }
        { "In " format.crossref.editor * }
	  if$
	      " \cite{" * crossref * "}" *
}

FUNCTION {in.book}
  { "In " output
      booktitle emphasize
	  volume empty$
	      {skip$}
          {", Vol.~" * volume *}
	    if$
	        edition empty$
		    {";" * output}
	        {"," * output
		       edition " ed.;\" * output}
		         if$
			     editor missing$
			         { skip$ }
			     { format.editors output }
			       if$
				   format.number.series output
				     new.sentence
				       publisher ":" * "publisher" output.check
				         address "," * "address" output.check
  }

FUNCTION {book.vol.ed}
{ "" volume empty$
      {skip$}
          {" Vol.~" * volume *}
	    if$
	        edition empty$
		    {";" * output}
	        {"," * output
		       edition " ed.;\" * output}
		         if$
}

FUNCTION {thesis.affilation}
{ school "," * "school" output.check
    institution empty$
          'skip$
	      {institution "," * output}
      if$
}

FUNCTION {article}
{ output.bibitem
    format.authors "author" output.check
        " " title * ", " * "title" output.check	
        crossref missing$
	    { journal field.or.null emphasize
	            format.yr.vol.pages
		            output
			        }
        { format.article.crossref output.nonnull
	        format.pages output
		      }
	  if$
	      note output
	        fin.entry
}

FUNCTION {submitted}
{ output.bibitem
    format.authors "author" output.check
        new.block
	  "submitted for publication in" output
	    note empty$
	      { journal emphasize "journal" output.check }
      { journal emphasizec "journal" output.check
	    note output }
	      if$
		  fin.entry
}

FUNCTION {inpress}
{ output.bibitem
    format.authors "author" output.check
        journal emphasizec "journal" output.check
	  "in press" output
	    note output
	      fin.entry
}

FUNCTION {book}
{ output.bibitem
    author empty$
          { format.editors "author and editor" output.check }
        { format.authors output.nonnull
	        crossref missing$
		     { "author and editor" editor either.or.check }
		   'skip$
		           if$
			       }
	  if$
	      booktitle empty$
	          { title emphasizesc output }
	      { booktitle emphasizesc output }
	        if$
		    crossref missing$
		        { format.bvolume output
			        format.number.series output
				        new.sentence
					      publisher ":" * "publisher" output.check
					            address "," * "address" output.check
						        }
		    { format.book.crossref output.nonnull }
		      if$
			  format.edition output
			    format.date "year" output.check
			      note output
			        fin.entry
}

FUNCTION {booklet}
{ output.bibitem
    format.authors "author" output.check
        " ``" title * "'',  " * "title" output.check
	  howpublished empty$
	      { skip$ }
        { howpublished "," * output}
	  if$
	      address empty$
	          { skip$ }
	      { address "," * output}
	        if$
		    format.date "year" output.check
		      note output
		        fin.entry
}

FUNCTION {inbook}
{ output.bibitem
    format.authors "author" output.check
        " " title * ".  " * "title" output.check
	  in.book
	    chapter empty$
	        {format.date "year" output.check}
        {year empty$
	        { "empty year field in " * cite$ * warning$ }
		      { year "; " * output
			        format.chapter.pages output }
				    if$
				          }
	  if$
	      note output
	        fin.entry
}

FUNCTION {incollection}
{ output.bibitem
    format.authors "author" output.check
        " " title * ".  " * "title" output.check
	  in.book
	    format.date "year" output.check
	      note output
	        fin.entry
}

FUNCTION {inproceedings}  { inbook }

FUNCTION {conference} { inbook }

FUNCTION {manual}
{ output.bibitem
    author empty$
          { organization empty$
	       'skip$
		    { organization ", " * address * output
		         }
	             if$
		           }
        { format.authors output.nonnull }
	  if$
	      "``" title * "'', " * "title" output.check
	        author empty$
		    { organization empty$
		         { address empty$
			               'skip$
					             { address "," * output }
				                 if$
						      }
			    'skip$
			            if$
				        }
	      { organization empty$
		        'skip$
			          { organization "," * output
				              address empty$
						              'skip$
							                    { address "," * output }
					                  if$
							              }
			     if$
			           }
	        if$
		    format.edition output
		      format.date output
		        note output
			  fin.entry
}

FUNCTION {mastersthesis}
{ output.bibitem
    format.authors "author" output.check
        "``" title * "'', " * "title" output.check
	  "Master's thesis, " output
	    thesis.affilation
	      format.date "year" output.check
	        note output
		  fin.entry
}

FUNCTION {misc}
{ output.bibitem
    format.authors output
        title empty$
	    'skip$
	        { "``" title * "'', " * output }
      if$
	  howpublished empty$
	      'skip$
	          {howpublished "," * output }
        if$
	    format.date output
	      note output
	        fin.entry
		  empty.misc.check
}

FUNCTION {phdthesis}
{ output.bibitem
    format.authors "author" output.check
        title emphasizec "title" output.check
	  format.thesis output.nonnull
	    thesis.affilation
	      format.date "year" output.check
	        note output
		  fin.entry
}

FUNCTION {proceedings}
{ output.bibitem
    editor empty$
          { organization output }
        { format.editors output.nonnull }
	  if$
	      booktitle empty$
	          { title emphasizec output }
	      { booktitle emphasizec output }
	        if$
		    book.vol.ed
		      format.number.series output
		        publisher empty$
			    { address empty$
			              { format.date "year" output.check  }
				              { address "," * output
						          format.date "year" output.check  }
							        if$
								        }
		    { address empty$
		              { publisher "," *  output
				          format.date "year" output.check  }
					          { publisher ": " * address * "," * output
						              format.date "year" output.check  }
							            if$
								            }
		      if$
			  note output
			    fin.entry
}

FUNCTION {techreport}
{ output.bibitem
    format.authors "author" output.check
        "``" title * "'', " * "title" output.check
	  type empty$
	      { "Technical Report" }
        'type
	    if$
	      number empty$
	           { ", " * }
	     { " " * number * "," * }
	       if$
		   output.nonnull
		     institution "," * "institution" output.check
		       address empty$
		           { skip$ }
	           { address "," * output }
		     if$
		         format.date "year" output.check
			   note output
			     fin.entry
}

FUNCTION {unpublished}
{ output.bibitem
    format.authors "author" output.check
        note empty$
	    { "Unpublished data" output }
        { year empty$
	         { note output }
		        { note ", " * year * output }
			      if$ }
			        if$
				    fin.entry
}

FUNCTION {remark}
{ output.bibitem
    note "note" output.check
        fin.entry
}

FUNCTION {default.type} { misc }

MACRO {jan} {"Jan."}
MACRO {feb} {"Feb."}
MACRO {mar} {"Mar."}
MACRO {apr} {"Apr."}
MACRO {may} {"May"}
MACRO {jun} {"June"}
MACRO {jul} {"July"}
MACRO {aug} {"Aug."}
MACRO {sep} {"Sept."}
MACRO {oct} {"Oct."}
MACRO {nov} {"Nov."}
MACRO {dec} {"Dec."}

MACRO {acchemr} {"Acc.\ Chem.\ Res."}
MACRO {aacsa} {"Adv.\ {ACS} Abstr."}
MACRO {anchem} {"Anal.\ Chem."}
MACRO {bioch} {"Biochemistry"}
MACRO {bicoc} {"Bioconj.\ Chem."}  % ***
MACRO {bitech} {"Biotechnol.\ Progr."}  % ***
MACRO {chemeng} {"Chem.\ Eng.\ News"}
MACRO {chs} {"Chem.\ Health Safety"} % ***
MACRO {crt} {"Chem.\ Res.\ Tox."} % ***
MACRO {chemrev} {"Chem.\ Rev."} % ***
MACRO {cmat} {"Chem.\ Mat."} % ***
MACRO {chemtech} {"{CHEMTECH}"} % ***
MACRO {enfu} {"Energy \&{} Fuels"} % ***
MACRO {envst} {"Environ.\ Sci.\ Technol."}
MACRO {iecf} {"Ind.\ Eng.\ Chem.\ Fundam."}
MACRO {iecpdd} {"Ind.\ Eng.\ Chem.\ Proc.\ Des.\ Dev."}
MACRO {iecprd} {"Ind.\ Eng.\ Chem.\ Prod.\ Res.\ Dev."}
MACRO {iecr} {"Ind.\ Eng.\ Chem.\ Res."} % ***
MACRO {inor} {"Inorg.\ Chem."}
MACRO {jafc} {"J.~Agric.\ Food Chem."}
MACRO {jacs} {"J.~Am.\ Chem.\ Soc."}
MACRO {jced} {"J.~Chem.\ Eng.\ Data"}
MACRO {jcics} {"J.~Chem.\ Inf.\ Comput.\ Sci."}
MACRO {jmc} {"J.~Med.\ Chem."}
MACRO {joc} {"J.~Org.\ Chem."}
MACRO {jps} {"J.~Pharm.\ Sci."}
MACRO {jpcrd} {"J.~Phys.\ Chem.\ Ref.\ Data"} % ***
MACRO {jpc} {"J.~Phys.\ Chem."}
MACRO {jpca} {"J.~Phys.\ Chem.~A"}
MACRO {jpcb} {"J.~Phys.\ Chem.~B"}
MACRO {lang} {"Langmuir"}
MACRO {macro} {"Macromolecules"}
MACRO {orgmet} {"Organometallics"}

MACRO {jft} {"J.~Chem.\ Soc., Faraday Trans."}
MACRO {jft1} {"J.~Chem.\ Soc., Faraday Trans. 1"}
MACRO {jft2} {"J.~Chem.\ Soc., Faraday Trans. 2"}
MACRO {tfs} {"Trans.\ Faraday Soc."}
MACRO {jcis} {"J.~Colloid Interface Sci."}
MACRO {acis} {"Adv.~Colloid Interface Sci."}
MACRO {cs} {"Colloids Surf."}
MACRO {csa} {"Colloids Surf.\ A:\ Physicochem.\ Eng.\ Apects"}
MACRO {csb} {"Colloids Surf.\ B:\ Biointerfaces"}
MACRO {pcps} {"Progr.\ Colloid Polym.\ Sci."}
MACRO {jmr} {"J.~Magn.\ Reson."}
MACRO {jmra} {"J.~Magn.\ Reson.\ A"}
MACRO {jmrb} {"J.~Magn.\ Reson.\ B"}
MACRO {sci} {"Science"}
MACRO {nat} {"Nature (London)"}
MACRO {jcch} {"J.~Comput.\ Chem."}
MACRO {cca} {"Croat.\ Chem.\ Acta"}

MACRO {poly} {"Polymer"}
MACRO {ajp} {"Am.\ J.\ Phys."}
MACRO {rsi} {"Rev.\ Sci.\ Instrum."}
MACRO {jcp} {"J.~Chem.\ Phys."}
MACRO {cpl} {"Chem.\ Phys.\ Lett."}
MACRO {molph} {"Mol.\ Phys."}
MACRO {pac} {"Pure Appl.\ Chem."}
MACRO {jbc} {"J.~Biol.\ Chem."}
MACRO {tl} {"Tetrahedron Lett.\"}
MACRO {psisoe} {"Proc.\ SPIE-Int.\ Soc.\ Opt.\ Eng."}
MACRO {prb} {"Phys.\ Rev.\ B:\ Condens.\ Matter"}
MACRO {jap} {"J.~Appl.\ Phys."}
MACRO {pnac} {"Proc.\ Natl.\ Acad.\ Sci.\ U.~S.~A. "}
MACRO {bba} {"Biochim.\ Biophys.\ Acta"}
MACRO {nar} {"Nucleic.\ Acid Res."}

READ

STRINGS { longest.label }

INTEGERS { number.label longest.label.width }

FUNCTION {initialize.longest.label}
{ "" 'longest.label :=
    #1 'number.label :=
    #0 'longest.label.width :=
}

FUNCTION {longest.label.pass}
{ number.label int.to.str$ 'label :=
    number.label #1 + 'number.label :=
        label width$ longest.label.width >
	    { label 'longest.label :=
	            label width$ 'longest.label.width :=
		          }
        'skip$
	    if$
}

EXECUTE {initialize.longest.label}

ITERATE {longest.label.pass}

FUNCTION {begin.bib}
{ preamble$ empty$
      'skip$
	    { preamble$ write$ newline$ }
        if$
	    "\providecommand{\refin}[1]{\\ \textbf{Referenced in:} #1}"
	      write$ newline$
	        "\begin{thebibliography}{"  longest.label  * "}" * write$ newline$
}

EXECUTE {begin.bib}

EXECUTE {init.state.consts}

ITERATE {call.type$}

FUNCTION {end.bib}
{ newline$
    "\end{thebibliography}" write$ newline$
}

EXECUTE {end.bib}
%% 
%%
%% End of file `achemso.bst'.
