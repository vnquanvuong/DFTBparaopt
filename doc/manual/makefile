tag=`date +%F`
inp=manual
out=${inp}_${tag}

all: 
	pdflatex   ${inp} 
	bibtex     ${inp} 
	makeindex  ${inp} 
	pdflatex   ${inp} 
	pdflatex   ${inp} 
	cp         ${inp}.pdf pdf/${out}.pdf

clean: 
	rm -rf *.pdf *.ilg *.aux *.brf *.ind *.log *.toc *.bbl *.idx *.out *blg *.tdo 
