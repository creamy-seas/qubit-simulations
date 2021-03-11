OSFLAG :=
ifeq ($(OS),Windows_NT)
	OSFLAG += WIN32
else
	OSFLAG += UNIX
endif

jupyter:
	@docker-compose up --build

# ilya-testing: $(OSFLAG)
	# @echo "Do not use - this is WIP for easier deployment"

UNIX:
	@echo "⚛  Assuming UNIX-like environment"
	@echo "⚙ Building system - this may take 5 minutes"
	@docker build -t qubit-simulations .
	@echo docker run --publish 8888:8888\
						--name qubit-simulations-container\
						--rm\
						-i\
						-t\
						--volume "$pwd":/home qubit-simulations

WIN:
	@echo "⚛  Assuming WINDOWS-like environment"
	@echo "⚙ Building system - this may take 5 minutes"
	@docker build -t qubit-simulations .
	@docker run --publish 8888:8888\
						--name qubit-simulations-container\
						--rm\
						-i\
						-t\
						--volume ${PWD}:/qutip-simulator\
						qubit-simulations
