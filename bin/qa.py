#!/usr/bin/python
import re
import sys
import hashlib
import subprocess as sp
import numpy as np

debug = False

typ1 = np.dtype([('name', 'a50'), ('beg', 'i'), ('end', 'i'), ('type', 'a4')])

# starts with spaces or spaces and one of { 'end', 'pure', ... }
# if function it can have a type next goes subroutine or function or type
test_for_routines = re.compile('''
      ^\s{0,12}(|end|pure|elemental|recursive|real|logical|integer)\s
      (|pure|elemental|recursive|real|logical|integer)(|\s)
      (subroutine|function|type(,|\s))
   ''', re.VERBOSE)
# starts with spaces or spaces and one of { 'end', 'pure', ... }
# next goes subroutine or function or type
test_for_interfaces = re.compile('''
      ^\s{0,12}(|end|abstract)\s
      interface
   ''', re.VERBOSE)
#test_for_routines  = re.compile('''
#       ^(?!\s{0,9}!).*(subroutine|function|type(,|\s::))
#   ''',re.VERBOSE)
module_body = re.compile(
    '''^\s{0,3}(module|contains|program)''', re.VERBOSE)
just_end = re.compile('''^\s{0,9}end''', re.IGNORECASE)

have_implicit = re.compile('''implicit\snone''', re.IGNORECASE)
have_privpub = re.compile('''^\s{0,9}(public|private)''', re.VERBOSE)
have_pub = re.compile('''^\s{0,9}public''', re.VERBOSE)
have_priv = re.compile('''^\s{0,9}private\s::''', re.VERBOSE)
remove_warn = re.compile('''(?!.*QA_WARN .+)''', re.VERBOSE)
have_global_public = re.compile('''^\s{0,9}public(?!.*::)''', re.VERBOSE)
depr_syntax_1 = re.compile('''
        ^\s{1,12}(?i)(?:real(?:\s|,)|integer(?:\s|,)|logical
        (?:\s|,|\()|character(?:\s|,))(?!.*::)
    ''', re.IGNORECASE)
depr_syntax_2 = re.compile(
    '''^\s{1,12}(?i)use[\s](?!.*only)''', re.IGNORECASE)
depr_syntax_3 = re.compile(
    '''^\s{1,12}(?i)character(?![(])''', re.IGNORECASE)
is_function = re.compile('''(?i)\sfunction\s''', re.IGNORECASE)
not_function = re.compile('''(?!.*function)''', re.IGNORECASE)
tab_char = re.compile('\t')
has_use = re.compile("^\s{1,12}(?i)use\s")
have_label = re.compile('^[0-9]', re.VERBOSE)
crude_write = re.compile("write *\( *\*", re.IGNORECASE)
magic_integer = re.compile("\(len=[1-9]", re.IGNORECASE)
continuation = re.compile('&$', re.VERBOSE)
implicit_save = re.compile('''
    (?i)(?:real(?:\s|,)|integer(?:\s|,)|logical(?:\s|,|\()|character(?:\s|,)).*
    ::.*=(|\s)[0-9]
''', re.IGNORECASE)
not_param_nor_save = re.compile("(?!.*(parameter|save))", re.IGNORECASE)

nasty_spaces = [
    re.compile("end\s{1,}do", re.IGNORECASE), "enddo",
    re.compile("end\s{1,}if", re.IGNORECASE), "endif",
    re.compile("end\s{1,}while", re.IGNORECASE), "endwhile",
    re.compile("end\s{1,}where", re.IGNORECASE), "endwhere",
    re.compile("only\s{1,}:", re.IGNORECASE), "only:",
    re.compile("if(|\s{2,})\(", re.IGNORECASE), "if (",
    re.compile("where(|\s{2,})\(", re.IGNORECASE), "where (",
    re.compile("while(|\s{2,})\(", re.IGNORECASE), "while (",
    re.compile("forall(|\s{2,})\(", re.IGNORECASE), "forall (",
    re.compile("\scase(|\s{2,})\(", re.IGNORECASE), " case ("
]


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''

b = bcolors()


def remove_binaries(files):
    list = []
    for file in files:
        checkFile = sp.Popen('file -bi ' + file, stdout=sp.PIPE, shell=True,
                             executable="/bin/bash")
        if not checkFile.communicate()[0].startswith('text'):
            print b.WARNING + "QA:  " + b.ENDC + file + \
                " is not a text file. I will not test it."
        else:
            list.append(file)
    return list


def select_sources(files):
    test = re.compile("F90$", re.IGNORECASE)
    return filter(test.search, files)


def wtf(lines, line, rname, fname):
    if(type(lines) == np.ndarray):
        linenum = line_num(lines, line)
    else:
        linenum = lines

    line = line.split("!")[0]    # Strip comments
    if(rname == ''):
        return " [%s]@L%i => %s" % (fname, linenum, line.strip())
    else:
        return " [%s:%s]@L%i => %s" % (fname, rname, linenum, line.strip())


def line_num(lines, line):
    return np.where(lines == line)[0][0]


def give_warn(s):
    return b.WARNING + s + b.ENDC


def give_err(s):
    return b.FAIL + s + b.ENDC


def parse_f90file(lines, fname, store):
    if (debug):
        print "[parse_f90file] fname = ", fname
    subs_array = np.zeros((0,), dtype=typ1)

    subs = filter(test_for_routines.search, lines)
    subs_names = []
    subs_types = []
    for f in subs:
        if (just_end.match(f)):
            word = f.strip().split(' ')
            subs_types.insert(0, word[1])
            subs_names.append(word[2])
    for f in subs_names:
        cur_sub = filter(re.compile(f).search, subs)
        if (len(cur_sub) > 2):
            if (debug):
                print "[parse_f90file] f, cur_sub = ", f, cur_sub
            for index in range(0, len(cur_sub)):
                if just_end.match(cur_sub[index]) and \
                    cur_sub[index][-len(f):] == f: break
        else:
            index = 1
        obj = (f, line_num(lines, cur_sub[index - 1]), line_num(
            lines, cur_sub[index]), subs_types.pop())
        subs_array = np.append(subs_array, np.array([obj], dtype=typ1))

    if (debug):
        print "[parse_f90file] subs = ", subs
        print "[parse_f90file] subs_names = ", subs_names

    mod = filter(module_body.match, lines)
    if (len(mod) > 1):
        obj = (mod[0].strip().split(" ")[1],
               line_num(lines, mod[0]),
               line_num(lines, mod[1]),
               mod[0].strip().split(" ")[0][0:3]
               )
        subs_array = np.append(subs_array, np.array([obj], dtype=typ1))
    elif (len(mod) == 1):
        obj = (mod[0].strip(
        ).split(" ")[1], line_num(lines, mod[0]), len(lines), 'mod')
        subs_array = np.append(subs_array, np.array([obj], dtype=typ1))
    else:
        store.append(
            give_warn("QA:  ") + "[%s] => module body not found!" % fname)
    return subs_array


def qa_check_id(store, fname):
    client = pysvn.Client()
    entry = client.info('.')

    try:
        for f in client.proplist(fname):
            pname, props = f
            fail = False
            if 'svn:keywords' in props:
                if 'Id' not in props['svn:keywords']:
                    fail = True
            else:
                fail = True
            if fail:
                store.append(
                    give_err("QA:  ") + "%s lacks svn:keywords Id" % (pname))
    except:
        return 0


def qa_checks(files, options):
    print b.OKBLUE + \
        '"I am the purifier, the light that clears all shadows."' + \
        ' - seal of cleansing inscription' + b.ENDC
    runpath = sys.argv[0].split("qa.py")[0]
    files = remove_binaries(files)
    # ToDo: check files other than F90
    f90files = select_sources(files)
    warns = []
    errors = []
    for f in f90files:
#        qa_check_id(errors, f)
        pfile = []
        lines = open(f, 'r').readlines()
        for line in lines:
            # things done in "in-place"
            line = line.rstrip()    # that removes trailing spaces
            for i in range(0, len(nasty_spaces), 2):
                line = re.sub(nasty_spaces[i], nasty_spaces[
                              i + 1], line)   # remove nasty spaces
            pfile.append(line)

        if lines != [line + '\n' for line in pfile]:
            fp = open(f, 'w')
            for line in pfile:
                fp.write(line + '\n')
            fp.close()

        #f = f.split('/')[-1]
        # checks for f90 file as whole
        qa_nonconforming_tabs(np.array(pfile), '', errors, f)
        qa_labels(np.array(pfile), '', errors, f)
        qa_crude_write(np.array(pfile), '', warns, f)
        qa_magic_integers(np.array(pfile), '', warns, f)
        # checks that require parsing f90 files
        clean_ind = []
        pfile = np.array(pfile)
        # remove interfaces as we currently don't handle them well
        interfaces = [line_num(
            pfile, i) for i in filter(test_for_interfaces.search, pfile)]
        while len(interfaces) > 0:
            if (debug):
                print "Removed interface"
            pfile = np.delete(pfile, np.s_[interfaces[0]:interfaces[1] + 1], 0)
            interfaces = [line_num(
                pfile, i) for i in filter(test_for_interfaces.search, pfile)]

        for obj in parse_f90file(pfile, f, warns):
            if (debug):
                print '[qa_checks] obj =', obj
            part = pfile[obj['beg']:obj['end']]
#         if (debug):
#            for f in part: print f
            # False refs need to be done before removal of types in module body
            qa_false_refs(part, obj['name'], warns, f)
            if(obj['type'] == 'mod'):
                # module body is always last, remove lines that've been
                # already checked
                part = np.delete(part, np.array(clean_ind) - obj['beg'])
                qa_have_priv_pub(part, obj['name'], warns, f)
            else:
                clean_ind += range(obj['beg'], obj['end'] + 1)

            qa_depreciated_syntax(part, obj['name'], warns, f)
            if(obj['type'] != 'type'):
                qa_have_implicit(part, obj['name'], errors, f)
                qa_implicit_saves(part, obj['name'], errors, f)

    for warning in warns:
        print warning
    for error in errors:
        print error

    if (len(errors)):
        print give_err("%i error(s) detected! " % len(
            errors)) + "I will not let you commit unless you force me!!!"
        if options.force:
            print "Damn! You are determined to make low quality commit :-("
        else:
            exit()
    else:
        print b.OKGREEN + "Yay! No errors!!! " + b.ENDC

    if (len(warns)):
        s = give_warn("%i warnings detected. " % len(warns)) + \
            "Do you wish to proceed? (y/N) "
        if(raw_input(s) != 'y'):
            print "See you later!"
            exit()
    else:
        print b.OKGREEN + "No warnings detected. " + b.ENDC + \
            "If everyone were like you, I'd be out of business!"


def qa_have_priv_pub(lines, name, warns, fname):
    if(not filter(have_privpub.search, lines)):
        warns.append(give_warn("QA:  ") +
                     "module [%s:%s] lacks public/private keywords." %
                     (fname, give_err(name)))
    else:
        if(filter(remove_warn.match, filter(have_priv.search, lines))):
            warns.append(give_warn("QA:  ") +
                         "module [%s:%s] have selective private." %
                         (fname, give_err(name)))
        if(filter(remove_warn.match,
                  filter(have_global_public.search, lines))):
            warns.append(give_warn("QA:  ") +
                         "module [%s:%s] is completely public." %
                         (fname, give_err(name)))


def qa_crude_write(lines, rname, store, fname):
    warning = 0
    for f in filter(remove_warn.match, filter(crude_write.search, lines)):
        store.append(
            give_warn("!! crude write  ") + wtf(lines, f, rname, fname))


def qa_magic_integers(lines, rname, store, fname):
    for f in filter(remove_warn.match, filter(magic_integer.search, lines)):
        hits = np.where(lines == f)[0]
        if(len(hits) > 1):
            for i in hits:
                warn = give_warn("!! magic integer") + wtf(i, f, rname, fname)
                if(warn not in store):
                    store.append(warn)
        else:
            warn = give_warn("!! magic integer") + wtf(lines, f, rname, fname)
            if(warn not in store):
                store.append(warn)


def qa_nonconforming_tabs(lines, rname, store, fname):
    for f in filter(tab_char.search, lines):
        store.append(give_err("QA:  ") + "non conforming tab detected " +
                     wtf(lines, f, rname, fname))


def qa_labels(lines, rname, store, fname):
    for f in filter(have_label.search, lines):
        store.append(give_err("QA:  ") + "label detected              " +
                     wtf(lines, f, rname, fname))


def qa_depreciated_syntax(lines, rname, store, fname):
#    print b.OKGREEN + "QA: " + b.ENDC + "Checking for depreciated syntax"
    for f in filter(not_function.match, filter(depr_syntax_1.search, lines)):
        store.append(
            give_warn("!! lacking ::   ") + wtf(lines, f, rname, fname))
    for f in filter(remove_warn.match, filter(depr_syntax_2.search, lines)):
        store.append(
            give_warn("!! greedy use   ") + wtf(lines, f, rname, fname))
    for f in filter(depr_syntax_3.search, lines):
        store.append(
            give_warn("!! wrong syntax ") + wtf(lines, f, rname, fname))


def qa_have_implicit(lines, name, store, fname):
    if(not filter(have_implicit.search, lines)):
        store.append(give_err(
            "QA:  ") + "missing 'implicit none'      [%s:%s]" % (fname, name))


def remove_amp(lines, strip):
    buf = ''
    temp = []
    for line in lines:
        if(len(buf)):
            line = buf + line.lstrip()
            buf = ''
        if(continuation.search(line)):
            buf = re.sub('&', '', line.split("!")[0])
        else:
            if(strip):
                temp.append(line.split("!")[0])  # kills QA_WARN
            else:
                temp.append(line)
    return temp


def qa_false_refs(lines, name, store, fname):
    temp = remove_amp(filter(remove_warn.match, lines), True)
    uses = filter(has_use.search, temp)

    for item in uses:
        to_check = [f.strip() for f in item.split("only:")[1].split(',')]
        to_check = [re.sub('&', '', f).lstrip(
        ) for f in to_check]     # additional sanitization
        # remove operator keyword from import
        for ino, item in enumerate(to_check):
            try:
                new_item = re.search('operator\((.+?)\)', item).group(1)
            except AttributeError:
                new_item = item
            to_check[ino] = new_item
        for func in to_check:
            pattern = re.compile(func, re.IGNORECASE)
            # stupid but seems to work
            if(len(filter(pattern.search, temp)) < 2):
                store.append(give_warn("QA:  ") + "'" + func +
                             "' grabbed but not used in [%s:%s]" %
                             (fname, name))


def qa_implicit_saves(lines, name, store, fname):
#   print b.OKGREEN + "QA: " + b.ENDC + "Checking for implicit saves"
    impl = filter(not_param_nor_save.match, filter(implicit_save.search,
                  remove_amp(filter(remove_warn.match, lines), True)))
    if(len(impl)):
        store.append(give_err(
            "QA:  ") + "implicit saves detected in   [%s:%s]" % (fname, name))
    for line in impl:
        store.append(line.strip())

if __name__ == "__main__":
    from optparse import OptionParser
    usage = "usage: %prog [options] FILES"
    parser = OptionParser(usage=usage)
    parser.add_option("-v", "--verbose",
                      action="store_true", dest="debug", default=False,
                      help="make lots of noise [default]")
    parser.add_option("-f", "--force",
                      action="store_true", dest="force",
                      help="commit despite errors (It will be logged)")
    (options, args) = parser.parse_args()

    debug = options.debug
    if len(args) < 1:
        parser.error("incorrect number of arguments")
    qa_checks(args, options)
