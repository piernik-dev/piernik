#/usr/bin/env bash

_setup_completions()
{
  if [ "${#COMP_WORDS[@]}" -gt "2" ]; then
    return
    # ToDo: add setup options here
  else
    COMPREPLY=($( compgen -W "$( find problems/ -type d | sed 's/^problems\///' ) --help --last --problems --units" -- "${COMP_WORDS[1]}" ))
  fi
}

complete -F _setup_completions ./setup
