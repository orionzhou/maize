#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import errno
import os
import time
import os.path as op
import shutil
import signal
import sys
import logging
import fnmatch

from urllib.parse import urlencode
from subprocess import PIPE, run

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self

def str2bool(v):
    if not isinstance(v, str):
        raise ValueError("invalid literal for boolean: Not a string.")
    if v.lower() in ("yes", "true", "t", "1"):
        return True
    elif v.lower() in ('no', 'false', 'f', '0'):
        return False
    else:
        raise ValueError('invalid literal for boolean: "%s"' % v)

def get_abs_path(link_name):
    source = link_name
    if op.islink(source):
        source = os.readlink(source)
    else:
        source = op.basename(source)

    link_dir = op.dirname(link_name)
    source = op.normpath(op.join(link_dir, source))
    source = op.abspath(source)
    if source == link_name:
        return source
    else:
        return get_abs_path(source)

def splitall(path):
    allparts = []
    while True:
        path, p1 = op.split(path)
        if not p1:
            break
        allparts.append(p1)
    allparts = allparts[::-1]
    return allparts

def get_module_docstring(filepath):
    "Get module-level docstring of Python module at filepath, e.g. 'path/to/file.py'."
    co = compile(open(filepath).read(), filepath, 'exec')
    if co.co_consts and isinstance(co.co_consts[0], basestring):
        docstring = co.co_consts[0]
    else:
        docstring = None
    return docstring

def dmain(mainfile, type="action"):
    cwd = op.dirname(mainfile)
    pyscripts = [x for x in glob(op.join(cwd, "*", '__main__.py'))] \
        if type == "module" \
        else glob(op.join(cwd, "*.py"))
    actions = []
    for ps in sorted(pyscripts):
        action = op.basename(op.dirname(ps)) \
            if type == "module" \
            else op.basename(ps).replace(".py", "")
        if action[0] == "_":  # hidden namespace
            continue
        pd = get_module_docstring(ps)
        action_help = [x.rstrip(":.,\n") for x in pd.splitlines(True) \
                if len(x.strip()) > 10 and x[0] != '%'][0] \
                if pd else "no docstring found"
        actions.append((action, action_help))

    a = ActionDispatcher(actions)
    a.print_help()

def backup(filename):
    bakname = filename + ".bak"
    if op.exists(filename):
        logging.debug("Backup `{0}` to `{1}`".format(filename, bakname))
        sh("mv {0} {1}".format(filename, bakname))
    return bakname

def sh(cmd, grid=False, infile=None, outfile=None, errfile=None,
        append=False, background=False, threaded=None, log=True,
        grid_opts=None, silent=False):
    """
    simple wrapper for system calls
    """
    if not cmd:
        return 1
    if silent:
        outfile = errfile = "/dev/null"
    # if grid:
        # from maize.apps.grid import GridProcess
        # pr = GridProcess(cmd, infile=infile, outfile=outfile, errfile=errfile,
                         # threaded=threaded, grid_opts=grid_opts)
        # pr.start()
        # return pr.jobid
    #else:
    if infile:
        cat = "cat"
        if infile.endswith(".gz"):
            cat = "zcat"
        cmd = "{0} {1} |".format(cat, infile) + cmd
    if outfile and outfile != "stdout":
        if outfile.endswith(".gz"):
            cmd += " | gzip"
        tag = ">"
        if append:
            tag = ">>"
        cmd += " {0}{1}".format(tag, outfile)
    if errfile:
        if errfile == outfile:
            errfile = "&1"
        cmd += " 2>{0}".format(errfile)
    if background:
        cmd += " &"

    if log:
        logging.debug(cmd)
    return run(cmd, shell=True)

def is_exe(fpath):
    return op.isfile(fpath) and os.access(fpath, os.X_OK)

def which(program):
    """
    Emulates the unix which command.

    >>> which("cat")
    "/bin/cat"
    >>> which("nosuchprogram")
    """
    fpath, fname = op.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = op.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def glob(pathname, pattern=None):
    """
    Wraps around glob.glob(), but return a sorted list.
    """
    import glob as gl
    if pattern:
        pathname = op.join(pathname, pattern)
    return natsorted(gl.glob(pathname))

def iglob(pathname, patterns):
    """
    Allow multiple file formats. This is also recursive. For example:

    >>> iglob("apps", "*.py,*.pyc")
    """
    matches = []
    patterns = patterns.split(",") if "," in patterns else listify(patterns)
    for root, dirnames, filenames in os.walk(pathname):
        matching = []
        for pattern in patterns:
            matching.extend(fnmatch.filter(filenames, pattern))
        for filename in matching:
            matches.append(op.join(root, filename))
    return natsorted(matches)

def symlink(target, link_name):
    try:
        os.symlink(target, link_name)
    except OSError as e:
        if e.errno == errno.EEXIST:
            os.remove(link_name)
            os.symlink(target, link_name)

def mkdir(dirname, overwrite=False):
    """
    Wraps around os.mkdir(), but checks for existence first.
    """
    if op.isdir(dirname):
        if overwrite:
            shutil.rmtree(dirname)
            os.mkdir(dirname)
            logging.debug("Overwrite folder `{0}`.".format(dirname))
        else:
            return False  # Nothing is changed
    else:
        try:
            os.mkdir(dirname)
        except:
            os.makedirs(dirname)
        logging.debug("`{0}` not found. Creating new.".format(dirname))

    return True

def is_newer_file(a, b):
    """
    Check if the file a is newer than file b
    """
    if not (op.exists(a) and op.exists(b)):
        return False
    am = os.stat(a).st_mtime
    bm = os.stat(b).st_mtime
    return am > bm

def parse_multi_values(param):
    values = None
    if param:
        if op.isfile(param):
            values = list(set(x.strip() for x in open(param)))
        else:
            values = list(set(param.split(",")))
    return values

def listify(a):
    return a if (isinstance(a, list) or isinstance(a, tuple)) else [a]

def last_updated(a):
    """
    Check the time since file was last updated.
    """
    return time.time() - op.getmtime(a)

def need_update(a, b):
    """
    Check if file a is newer than file b and decide whether or not to update
    file b. Can generalize to two lists.
    """
    a = listify(a)
    b = listify(b)

    return any((not op.exists(x)) for x in b) or \
           all((os.stat(x).st_size == 0 for x in b)) or \
           any(is_newer_file(x, y) for x in a for y in b)

def get_today():
    """
    Returns the date in 2010-07-14 format
    """
    from datetime import date
    return str(date.today())

def ls_ftp(dir):
    from urlparse import urlparse
    from ftplib import FTP, error_perm
    o = urlparse(dir)

    ftp = FTP(o.netloc)
    ftp.login()
    ftp.cwd(o.path)

    files = []
    try:
        files = ftp.nlst()
    except error_perm as resp:
        if str(resp) == "550 No files found":
            print("no files in this directory")
        else:
            raise
    return files

def download(url, filename=None, debug=True, cookies=None):
    from urlparse import urlsplit
    from subprocess import CalledProcessError
    from jcvi.formats.base import FileShredder

    scheme, netloc, path, query, fragment = urlsplit(url)
    filename = filename or op.basename(path)
    filename = filename.strip()

    if not filename:
        filename = "index.html"

    if op.exists(filename):
        if debug:
            msg = "File `{0}` exists. Download skipped.".format(filename)
            logging.error(msg)
    else:
        from jcvi.utils.ez_setup import get_best_downloader

        downloader = get_best_downloader()
        try:
            downloader(url, filename, cookies=cookies)
        except (CalledProcessError, KeyboardInterrupt) as e:
            print >> sys.stderr, e
            FileShredder([filename])

    return filename

def getfilesize(filename, ratio=None):
    rawsize = op.getsize(filename)
    if not filename.endswith(".gz"):
        return rawsize

    import struct

    fo = open(filename, 'rb')
    fo.seek(-4, 2)
    r = fo.read()
    fo.close()
    size = struct.unpack('<I', r)[0]
    # This is only ISIZE, which is the UNCOMPRESSED modulo 2 ** 32
    if ratio is None:
        return size

    # Heuristic
    heuristicsize = rawsize / ratio
    while size < heuristicsize:
        size += 2 ** 32
    if size > 2 ** 32:
        logging.warn(\
            "Gzip file estimated uncompressed size: {0}.".format(size))

    return size

def debug():
    from maize.apps.console import magenta, yellow

    format = yellow("%(asctime)s [%(module)s]")
    format += magenta(" %(message)s")
    logging.basicConfig(level=logging.DEBUG,
            format=format,
            datefmt="%H:%M:%S")

debug()

def expand(args):
    """
    %prog expand */*

    Move files in subfolders into the current folder. Use --symlink to 
    create a link instead.
    """
    seen = set()
    for a in args.dirs:
        oa = a.replace("/", "_")
        if oa in seen:
            logging.debug("Name collision `{0}`, ignored.".format(oa))
            continue

        cmd = "cp -s" if args.symlink else "mv"
        cmd += " {0} {1}".format(a, oa)
        sh(cmd)
        seen.add(oa)

def mvfolder(args):
    dirw = args.dirw
    print("working in %s" % dirw)
    os.chdir(dirw)
    ary = os.listdir(dirw)
    for item in ary:
        print(item)
        if item == 'temp' or item.startswith("sna"):
            continue
        cmd = "cp -rf %s %s.bak" % (item, item)
        print("  " + cmd)
        os.system(cmd)
        cmd = "rm -rf %s" % item
        print("  " + cmd)
        os.system(cmd)
        cmd = "mv %s.bak %s" % (item, item)
        print("  " + cmd)
        os.system(cmd)

def fname():
    return sys._getframe().f_back.f_code.co_name

def get_times(filename):
    st = os.stat(filename)
    atime = st.st_atime
    mtime = st.st_mtime
    return (atime, mtime)

def timestamp(args):
    """
    %prog timestamp path > timestamp.info

    Record the timestamps for all files in the current folder.
    filename	atime	mtime

    This file can be used later to recover previous timestamps through touch().
    """
    path = args.dir
    for root, dirs, files in os.walk(path):
        for f in files:
            filename = op.join(root, f)
            atime, mtime = get_times(filename)
            print(filename, atime, mtime)

def touch(args):
    """
    %prog touch timestamp.info

    Recover timestamps for files in the current folder.
    CAUTION: you must execute this in the same directory as timestamp().
    """
    from time import ctime

    info = args.dir
    fp = open(info)
    for row in fp:
        path, atime, mtime = row.split()
        atime = float(atime)
        mtime = float(mtime)
        current_atime, current_mtime = get_times(path)

        # Check if the time has changed, with resolution up to 1 sec
        if int(atime) == int(current_atime) and \
           int(mtime) == int(current_mtime):
            continue

        times = [ctime(x) for x in (current_atime, current_mtime, atime, mtime)]
        msg = "{0} : ".format(path)
        msg += "({0}, {1}) => ({2}, {3})".format(*times)
        print >> sys.stderr, msg
        os.utime(path, (atime, mtime))

def snapshot(fp, p, fsize, counts=None):
    pos = int(p * fsize)
    print("==>> File `{0}`: {1} ({2}%)".format(fp.name, pos, int(p * 100)))
    fp.seek(pos)
    fp.next()
    for i, row in enumerate(fp):
        if counts and i > counts:
            break
        try:
            sys.stdout.write(row)
        except IOError:
            break

def less(args):
    """
    %prog less filename position | less

    Enhance the unix `less` command by seeking to a file location first. This is
    useful to browse big files. Position is relative 0.00 - 1.00, or bytenumber.

    $ %prog less myfile 0.1      # Go to 10% of the current file and streaming
    $ %prog less myfile 0.1,0.2  # Stream at several positions
    $ %prog less myfile 100      # Go to certain byte number and streaming
    $ %prog less myfile 100,200  # Stream at several positions
    $ %prog less myfile all      # Generate a snapshot every 10% (10%, 20%, ..)
    """
    from jcvi.formats.base import must_open

    filename, pos = args.fi, args.pos
    fsize = getfilesize(filename)

    if pos == "all":
        pos = [x / 10. for x in range(0, 10)]
    else:
        pos = [float(x) for x in pos.split(",")]

    if pos[0] > 1:
        pos = [x / fsize for x in pos]

    if len(pos) > 1:
        counts = 20
    else:
        counts = None

    fp = must_open(filename)
    for p in pos:
        snapshot(fp, p, fsize, counts=counts)

def send_email(fromaddr, toaddr, subject, message):
    """
    Send an email message
    """
    from smtplib import SMTP
    from email.mime.text import MIMEText

    SERVER = "localhost"
    _message = MIMEText(message)
    _message['Subject'] = subject
    _message['From'] = fromaddr
    _message['To'] = ", ".join(toaddr)

    server = SMTP(SERVER)
    server.sendmail(fromaddr, toaddr, _message.as_string())
    server.quit()

def is_valid_email(email):
    """
    RFC822 Email Address Regex
    --------------------------

    Originally written by Cal Henderson
    c.f. http://iamcal.com/publish/articles/php/parsing_email/

    Translated to Python by Tim Fletcher, with changes suggested by Dan Kubb.

    Licensed under a Creative Commons Attribution-ShareAlike 2.5 License
    http://creativecommons.org/licenses/by-sa/2.5/
    """
    import re

    qtext = '[^\\x0d\\x22\\x5c\\x80-\\xff]'
    dtext = '[^\\x0d\\x5b-\\x5d\\x80-\\xff]'
    atom = '[^\\x00-\\x20\\x22\\x28\\x29\\x2c\\x2e\\x3a-\\x3c\\x3e\\x40\\x5b-\\x5d\\x7f-\\xff]+'
    quoted_pair = '\\x5c[\\x00-\\x7f]'
    domain_literal = "\\x5b(?:%s|%s)*\\x5d" % (dtext, quoted_pair)

    quoted_string = "\\x22(?:%s|%s)*\\x22" % (qtext, quoted_pair)
    domain_ref = atom
    sub_domain = "(?:%s|%s)" % (domain_ref, domain_literal)
    word = "(?:%s|%s)" % (atom, quoted_string)
    domain = "%s(?:\\x2e%s)*" % (sub_domain, sub_domain)
    local_part = "%s(?:\\x2e%s)*" % (word, word)
    addr_spec = "%s\\x40%s" % (local_part, domain)

    email_address = re.compile('\A%s\Z' % addr_spec)
    if email_address.match(email):
        return True
    return False

def pid_exists(pid):
    """Check whether pid exists in the current process table."""
    if pid < 0:
        return False
    import errno
    try:
        os.kill(pid, 0)
    except OSError as e:
        return e.errno == errno.EPERM
    else:
        return True

class TimeoutExpired(Exception):
    pass

def _waitpid(pid, interval=None, timeout=None):
    """
    Wait for process with pid 'pid' to terminate and return its
    exit status code as an integer.

    If pid is not a children of os.getpid() (current process) just
    waits until the process disappears and return None.

    If pid does not exist at all return None immediately.

    Raise TimeoutExpired on timeout expired (if specified).

    Source: http://code.activestate.com/recipes/578022-wait-for-pid-and-check-for-pid-existance-posix
    """
    def check_timeout(delay):
        if timeout is not None:
            if time.time() >= stop_at:
                raise TimeoutExpired
        time.sleep(delay)
        return min(delay * 2, interval)

    if timeout is not None:
        waitcall = lambda: os.waitpid(pid, os.WNOHANG)
        stop_at = time.time() + timeout
    else:
        waitcall = lambda: os.waitpid(pid, 0)

    delay = 0.0001
    while 1:
        try:
            retpid, status = waitcall()
        except OSError as err:
            if err.errno == errno.EINTR:
                delay = check_timeout(delay)
                continue
            elif err.errno == errno.ECHILD:
                # This has two meanings:
                # - pid is not a child of os.getpid() in which case
                #   we keep polling until it's gone
                # - pid never existed in the first place
                # In both cases we'll eventually return None as we
                # can't determine its exit status code.
                while 1:
                    if pid_exists(pid):
                        delay = check_timeout(delay)
                    else:
                        return
            else:
                raise
        else:
            if retpid == 0:
                # WNOHANG was used, pid is still running
                delay = check_timeout(delay)
                continue

        # process exited due to a signal; return the integer of
        # that signal
        if os.WIFSIGNALED(status):
            return os.WTERMSIG(status)
        # process exited using exit(2) system call; return the
        # integer exit(2) system call has been called with
        elif os.WIFEXITED(status):
            return os.WEXITSTATUS(status)
        else:
            # should never happen
            raise RuntimeError("unknown process exit status")

def waitpid(args):
    """
    %prog waitpid PID ::: "./command_to_run param1 param2 ...."

    Given a PID, this script will wait for the PID to finish running and
    then perform a desired action (notify user and/or execute a new command)

    Specify "--notify=METHOD` to send the user a notification after waiting for PID
    Specify `--grid` option to send the new process to the grid after waiting for PID
    """
    import shlex
    from jcvi.utils.iter import flatten

    valid_notif_methods.extend(list(flatten(available_push_api.values())))

    p = OptionParser(waitpid.__doc__)
    p.add_option("--notify", default="email", choices=valid_notif_methods,
                 help="Specify type of notification to be sent after waiting")
    p.add_option("--interval", default=120, type="int",
                 help="Specify PID polling interval in seconds")
    p.add_option("--message",
                help="Specify notification message [default: %default]")
    p.set_email()
    p.set_grid()
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())

    if not opts.message:
        """
        If notification message not specified by user, just get
        the name of the running command and use it as the message
        """
        from subprocess import check_output

    sep = ":::"
    cmd = None
    if sep in args:
        sepidx = args.index(sep)
        cmd = " ".join(args[sepidx + 1:]).strip()
        args = args[:sepidx]

    pid = int(" ".join(args).strip())

    status = pid_exists(pid)
    if status:
        if opts.message:
            msg = opts.message
        else:
            get_origcmd = "ps -p {0} -o cmd h".format(pid)
            msg = check_output(shlex.split(get_origcmd)).strip()
        _waitpid(pid, interval=opts.interval)
    else:
        logging.debug("Process with PID {0} does not exist".format(pid))
        sys.exit()

    if opts.notify:
        notifycmd = ["[{0}] `{1}`".format(gethostname(), msg)]
        if opts.notify != "email":
            notifycmd.append("--method={0}".format("push"))
            notifycmd.append("--api={0}".format(opts.notify))
        else:
            notifycmd.append('--email={0}'.format(opts.email))
        notify(notifycmd)

    if cmd is not None:
        bg = False if opts.grid else True
        sh(cmd, grid=opts.grid, background=bg)

def inspect(object):
    """ A better dir() showing attributes and values
    """
    for k in dir(object):
        try:
            details = getattr(object, k)
        except Exception as e:
            details = e

        try:
            details = str(details)
        except Exception as e:
            details = e

        print >> sys.stderr, "{}: {}".format(k, details)

def sample_N(a, N):
    """ When size of N is >= size of a, random.sample() will emit an error:
    ValueError: sample larger than population

    This method handles such restrictions by repeatedly sampling when that
    happens.

    Examples:
    >>> sample_N([1, 2, 3], 2)
    >>> sample_N([1, 2, 3], 3)
    >>> sample_N([1, 2, 3], 4)
    """
    import random

    if N < len(a):
        return random.sample(a, N)

    return [random.choice(a) for x in range(N)]

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog = 'python -m maize.apps.base',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'basic support for running library as script'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("less", help = "enhanced version of the unix `less` command")
    sp1.add_argument('fi', help = 'input file')
    sp1.add_argument('pos', help = 'file position')
    sp1.set_defaults(func = less)

    sp1 = sp.add_parser("timestamp", help = "record timestamps for all files in the current folder")
    sp1.add_argument('dir', help = 'directory path')
    sp1.set_defaults(func = timestamp)

    sp1 = sp.add_parser("mvfolder", help='recursively move files/folders',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('dirw', nargs='?', default='/scratch.global/zhoux379', help='workding direcotry')
    sp1.set_defaults(func = mvfolder)

    sp1 = sp.add_parser("expand", help = "move files in subfolders into the current folder")
    sp1.add_argument('subdir', nargs = "+", help = 'sub-directory path')
    sp1.add_argument("--symlink", action="store_true", help="create symlink")
    sp1.set_defaults(func = expand)

    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()

