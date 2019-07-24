::latex --synctex=1 OneReport
::bibtex OneReport
::latex --synctex=1 OneReport
::latex --synctex=1 OneReport
::dvipdfm OneReport

pdflatex --synctex=1 OneReport >nul
::bibtex OneReport >nul
::pdflatex --synctex=1 OneReport >nul
@pdflatex --synctex=1 OneReport >nul
::dvipdfm OneReport
