import contextvars

log_file_var = contextvars.ContextVar("log_file", default = "gtotree-runlog.txt")
