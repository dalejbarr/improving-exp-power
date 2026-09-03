#!/bin/bash
#
# this script retrieved from: https://gist.github.com/padawanphysicist/d6299870de4ef8ad892f
# 
# I wrapped the code constructed in
#
# http://hacks-galore.org/aleix/blog/archives/2013/01/08/install-emacs-packages-from-command-line
#
# in a single bash script, so I would a single code snippet.
#

# Package to be installed
pkg_name=$1

# Elisp script is created as a temporary file, to be removed after installing 
# the package
elisp_script_name=$(mktemp /tmp/emacs-pkg-install-el.XXXXXX)
elisp_code="
;;
;; Install package from command line. Example:
;;
;;   $ emacs --batch --expr \"(define pkg-to-install 'smex)\" -l emacs-pkg-install.el
;;
(require 'package)
(setq package-user-dir (car package-directory-list)) ;; \"/opt/emacs.d/packages\"
(let* ((no-ssl (and (memq system-type '(windows-nt ms-dos))
                    (not (gnutls-available-p))))
       (proto (if no-ssl \"http\" \"https\")))
  (when no-ssl
    (warn \"\
	  Your version of Emacs does not support SSL connections,
which is unsafe because it allows man-in-the-middle attacks.
There are two things you can do about this warning:
1. Install an Emacs version that does support SSL and be safe.
2. Remove this warning from your init file so you won't see it again.\"))
  ;; Comment/uncomment these two lines to enable/disable MELPA and MELPA Stable as desired
  (add-to-list 'package-archives (cons \"melpa\" (concat proto \"://melpa.org/packages/\")) t)
  ;;(add-to-list 'package-archives (cons \"melpa-stable\" (concat proto \"://stable.melpa.org/packages/\")) t)
  (when (< emacs-major-version 24)
    ;; For important compatibility libraries like cl-lib
    (add-to-list 'package-archives (cons \"gnu\" (concat proto \"://elpa.gnu.org/packages/\")))))
;; (add-to-list 'package-archives
;;             '(\"marmalade\" . \"http://marmalade-repo.org/packages/\") t)
(package-initialize)
;; Fix HTTP1/1.1 problems
(setq url-http-attempt-keepalives nil)
(package-refresh-contents)
(package-install pkg-to-install)"

echo "$elisp_code" > $elisp_script_name

if [ $# -ne 1 ]
then
  echo "Usage: `basename $0` <package>"
  exit 1
fi

emacs --batch --eval "(defconst pkg-to-install '$pkg_name)" -l $elisp_script_name

# Remove tmp file
rm "$elisp_script_name"
