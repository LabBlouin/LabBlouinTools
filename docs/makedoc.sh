epydoc --html --name LabBlouinTools --url "http://LabBlouin.github.io/LabBlouinTools/" --graph all *.py -o /githubpages/LabBlouinTools/
epydoc --pdf --name LabBlouinTools -o ../docs/ *.py
find ../docs/ -type f -not -name api.pdf -not -name makedoc.sh -delete