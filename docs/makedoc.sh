epydoc --pdf --name LabBlouinTools -o ../docs/ *.py
find ../docs/ -type f -not -name api.pdf -not -name makedoc.sh -delete
cd .. ; sphinx-apidoc -F -e -H LabBlouinTools -A "Christian Blouin, Alex Safatli, Jose Sergio Hleap" -V $1 -R $1 -o docs labblouin ; cd labblouin ; rm *.pyc
cd ../docs/ ; make html ; cd -