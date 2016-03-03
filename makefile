################################################################
## Management tasks for the project

## This is a target to print some help message
hello:
	@echo "Hello"

help:
	@echo "Supported targets:"
	@echo "	github	open a connection to the github Web site"


github:
	@open -a firefox https://github.com/jvanheld/rna-seq_replicate_analysis

