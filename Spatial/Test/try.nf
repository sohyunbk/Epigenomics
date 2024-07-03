process myProcess {
    publishDir '/path/to/outputdir', mode: 'copy'  // or 'symlink', 'move'

    script:
    """
    # Your script commands
    """
}
