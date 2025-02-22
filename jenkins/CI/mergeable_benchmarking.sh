date
git remote add gawrysz https://github.com/gawrysz/piernik
git fetch gawrysz
git log --graph --branches --decorate --oneline -30
git checkout -b ___bm
git config --global user.email "Jenkins@no.email"
git config --global user.name "Jenkins test"
git merge remotes/gawrysz/benchmarking_mergeable | tee merge.txt
[ $( git ls-files -u | wc -l ) -le 0 ] || echo "Warning : Cannot merge cleanly with benchmarking_mergeable"

