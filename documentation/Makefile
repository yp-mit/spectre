USE_PDFLATEX=1


#CVTOPT_skinheads.gif := -density 300x300

#export Path := $(PATH);$(GSDIRB)\lib;$(GSDIRB)\bin;$(TEXDIRB)
#export PATH := $(Path)

figures_eps := $(wildcard plots/*.eps)

source_figs := $(wildcard figs/*.fig)
figures_all := $(figures_eps:.eps=.pdf) $(source_figs:.fig=.pdf)

### Conditionally add dependence on *.bbl file
ifndef NOBIBTEX
manual.pdf: manual.bbl
%.bbl: refs.bib %.tex
	if [ ! -e $*.aux ]; then make NOBIBTEX=1 $*.pdf; fi
	if ! fgrep '\bibdata' $*.aux; then touch $*.tex; make NOBIBTEX=1 $*.pdf; fi
	bibtex $*
endif

manual.pdf: manual.tex $(figures_all)

spellcheck:
	for fname in ${source_texs}; do \
		aspell --mode=tex --ignore 2 \
			--add-tex-command="eqref p" --add-tex-command="sqrt p" \
			--add-tex-command="log p" check $$fname; \
	done

figs/%.eps figs/%.pdf: figs/%.fig figs/%.tex
	fig2dev -L eps figs/$*.fig figs/$*_tmp.eps
	cd figs; latex $*.tex 
	cd figs; dvips $*.dvi -o $*.ps
	ps2eps -P -f figs/$*.ps
	@export BB1 BB2 BB3 BB4; egrep "^%%BoundingBox:" "figs/$*.eps" | while read DUMMY BB1 BB2 BB3 BB4; do \
		echo "======= Converting eps to pdf ========" ; \
		echo "    Source file     : figs/$*.eps " ;	\
		echo "    Destination file: figs/$*.pdf " ;	\
		echo "    Bounding box    : $$BB1 $$BB2 $$BB3 $$BB4"; \
		echo "======================================" ; \
		gs -dAutoRotatePages=/None -sDEVICE=pdfwrite -o "figs/$*.pdf" -c \
			"<< /PageSize [$$BB3 $$BB1 sub $$BB4 $$BB2 sub ] >> setpagedevice" \
			0 $$BB1 sub 0 $$BB2 sub translate -f "figs/$*.eps";  \
		break; \
	done;
#	rm -f figs/$*_tmp.eps figs/$*.dvi figs/$*.ps figs/$*.eps


%.fig: %.eps
	pstoedit -dis -f fig $< >$@
#### Special rules to create pdf images
%.pdf: %.eps
	@export BB1 BB2 BB3 BB4; egrep -a "^%%BoundingBox:" "$<" | while read DUMMY BB1 BB2 BB3 BB4; do \
		echo "======= Converting eps to pdf ========" ; \
		echo "    Source file     : $< " ;	\
		echo "    Destination file: $@ " ;	\
		echo "    Bounding box    : $$BB1 $$BB2 $$BB3 $$BB4"; \
		echo "======================================" ; \
		gs -dAutoRotatePages=/None -sDEVICE=pdfwrite -o "$@" -c \
			"<< /PageSize [$$BB3 $$BB1 sub $$BB4 $$BB2 sub ] >> setpagedevice" \
			0 $$BB1 sub 0 $$BB2 sub translate -f "$<";  \
	done;

#pages: $(object:.pdf=.ps)
#	$(GSDIRB)\bin\gswin32c -dSAFER -dBATCH -dNOPAUSE -r150 -sDEVICE=pnggray -dTextAlphaBits=4 -sOutputFile=pages-%02d.png $(object:.pdf=.ps)

ifdef USE_PDFLATEX
%.pdf: %.tex
	pdflatex -src-specials $<
	gs -q -dSAFER -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -dPDFSETTINGS=/printer \
		-dCompatibilityLevel=1.4 -sPAPERSIZE=letter -sOutputFile=$*_out.pdf -f $@
	mv $*_out.pdf $@
endif # USE_PDFLATEX #

%.pdf: %.ps
#
#  Options to consider in ps2pdf:
#
#  -dPDFSETTINGS=configuration 
#	Presets the "distiller parameters" to one of four predefined settings: 
#	/screen selects low-resolution output similar to the Acrobat Distiller "Screen Optimized" 
#		setting. 
#	/printer selects output similar to the Acrobat Distiller "Print Optimized" setting. 
#		NOTE: this setting also embeds ALL fonts as required by EDAS (e.g.)
#	/prepress selects output similar to Acrobat Distiller "Prepress Optimized" setting. 
#	/default selects output intended to be useful across a wide variety of uses, possibly 
#		at the expense of a larger output file.
#
#  -r<resolution>
#  	Consider trying -r8000 if encounter any problems
#
	gs -q -dSAFER -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -dPDFSETTINGS=/printer \
	-dCompatibilityLevel=1.4 -sPAPERSIZE=letter -sOutputFile=$@ -f $<

%.ps: %.dvi
#	DVIPS options:
# 	-j0 	turn off "partial fonts downloading (pfd)". PFD means that dvips embeds only
#  		glyphs that were actually used. Older dvips sometimes hangs or bugs w/o -j0
#  		
#  	-Ppdf	sets DPI to 8000 (so some warnings under tetex-2.0, I removed this under linux)
#
#  	-Pcm-super	use CM-SUPER fonts (in linux used by default anyway)
	
	dvips  $< -o $@
	
%.dvi: %.tex
	latex -src-specials $<


clean:
	rm -f $(objects_ps) $(objects_pdf) $(objects_dvi) \
	$(foreach obj, $(objects_pdf), $(obj:.pdf=.toc)) \
	$(foreach obj, $(objects_pdf), $(obj:.pdf=.aux)) \
	$(foreach obj, $(objects_pdf), $(obj:.pdf=.lof)) \
	$(foreach obj, $(objects_pdf), $(obj:.pdf=.log)) \
	$(foreach obj, $(objects_pdf), $(obj:.pdf=.nav)) \
	$(foreach obj, $(objects_pdf), $(obj:.pdf=.out)) \
	$(foreach obj, $(objects_pdf), $(obj:.pdf=.snm)) \
	$(foreach fig, $(source_figs), $(fig:.fig=.eps)) \
	$(foreach fig, $(source_figs), $(fig:.fig=.ps)) \
	$(foreach fig, $(source_figs), $(fig:.fig=.pdf)) \
	$(foreach fig, $(source_figs), $(fig:.fig=)_tmp.eps) \
	$(foreach fig, $(source_figs), $(fig:.fig=.aux)) \
	$(foreach fig, $(source_figs), $(fig:.fig=.log)) \
	$(foreach fig, $(source_figs), $(fig:.fig=.dvi)) \
	missfont.log 
	

%.eps: %.gif
	convert $< $(CVTOPT_$<) $@

%.eps: %.svg
	convert $< $(CVTOPT_$<) $@

%.eps: %.png
	convert $< $(CVTOPT_$<) $@
	
%.eps: %.bmp
	convert $< $(CVTOPT_$<) $@

%.eps: %.tif
	convert $< $(CVTOPT_$<) $@

%.eps: %.jpg
	convert $< $(CVTOPT_$<) $@
	
%.pdf: %.gif
	convert $< $(CVTOPT_$<) $@
	
%.pdf: %.bmp
	convert $< $(CVTOPT_$<) $@
