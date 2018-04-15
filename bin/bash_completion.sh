#/usr/bin/env bash

# Call 'source bin/bash_completion.sh' to get bash autocompletion of problem names and setup options.
# You can add it to your bash startup files too.

_piernik_completions()
{
  COMPREPLY=()
  if [ "${#COMP_WORDS[@]}" -gt "2" ]; then
    COMPREPLY=($( compgen -W "-v --verbose -q --laconic --debug \
    -n --nocompile --copy -l --linkexe -p --param= -d --define= --f90flags= \
    -c --compiler= -o --obj=" -- "${COMP_WORDS[COMP_CWORD]}" ))
    [[ $COMPREPLY == *= ]] && compopt -o nospace
  else
    COMPREPLY=($( compgen -W "$( find problems/ -type d | sed 's/^problems\///' ) \
    -h --help --last --problems -u --units" -- "${COMP_WORDS[1]}" ))
  fi
}

for i in ./setup piernik_setup.py ; do
  complete -F _piernik_completions $i
done
