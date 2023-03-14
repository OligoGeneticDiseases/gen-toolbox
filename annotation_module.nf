// Define the 'findType' process that takes in a command and arguments
process findType {
    input:
    val command
    val args
    // The script that will be executed for this process, using the input command and arguments
    script:
    """
    python main.py ${command} ${args}
    """
}

// Define the 'readVCFS' process that takes in a command and arguments
process readVCFS {
    input:
    val command
    val args
    // The script that will be executed for this process, using the input command and arguments
    script:
    """
    python main.py ${command} ${args}
    """
}

// Define the 'loadDB' process that takes in a command and arguments
process loadDB {
    input:
    val command
    val args
    // The script that will be executed for this process, using the input command and arguments
    script:
    """
    python main.py ${command} ${args}
    """
}

// Define the main workflow
workflow {
    // Check if the 'command' parameter was passed to the workflow
    if (!params.command) {
        println "Command not found"
        exit 1
    }
    // If the command is 'findtype', execute the 'findType' process with the input command and arguments
    if (params.command == "findtype") {
        findType(command=params.command, args=params.args)
    }
    // If the command is 'readvcfs', execute the 'readVCFS' process with the input command and arguments
    else if (params.command == "readvcfs") {
        readVCFS(command=params.command, args=params.args)
    }
    // If the command is 'loaddb', execute the 'loadDB' process with the input command and arguments
    else if (params.command == "loaddb") {
        loadDB(command=params.command, args=params.args)
    }
}