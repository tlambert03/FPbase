lint:
	flake8 --ignore=W503,I100,I101,I201,I202 favit fpbase fpseq proteins references

format:
	black --exclude "/migrations/" favit fpbase fpseq proteins references