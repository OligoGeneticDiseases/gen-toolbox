// Define the input and output channels for your module
input:
  file source_dir
  file dest_dir
  string file_type

output:
  file "output.txt"

// Define the command that will run your Python code
command = """
python main.py findtype -s ${source_dir} -d ${dest_dir} -t ${file_type}
"""

// Invoke the command using the Nextflow process directive
process runCommand {
  input:
    file source_dir
    file dest_dir
    string file_type

  output:
    file "output.txt"

  script:
    // Load the utility functions and command creator/handler classes
    include {utils.py}
    include {CommandFactory.py}
    include {CommandHandler.py}

    // Invoke the command using the command handler
    command = CommandFactory.createFindTypeCommand(source_dir, dest_dir, file_type)
    output = CommandHandler.invokeCommand(command)

    // Write the output to a file
    writeFile("output.txt", output)
}