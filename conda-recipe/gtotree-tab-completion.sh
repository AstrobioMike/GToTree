#!/usr/bin/env bash

### tab completion for the `gtt` command via argcomplete ###

# only activate in interactive shells
case "$-" in
    *i*) ;;
    *) return 0 2>/dev/null || exit 0 ;;
esac

if command -v register-python-argcomplete >/dev/null 2>&1; then
    eval "$(register-python-argcomplete gtt)"
fi
