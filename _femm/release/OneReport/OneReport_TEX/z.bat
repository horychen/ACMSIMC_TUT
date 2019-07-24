latex --synctex=1 OneReport
::bibtex OneReport >nul
latex --synctex=1 OneReport >nul
::latex --synctex=1 OneReport >nul
dvipdfm OneReport