# NewBkgMethodExtendedSource
New background method for extended source analyses

  Some simple git commands:
  
  How to create a remote branch:
  
          First, you create your branch locally:
  
          git checkout -b <branch-name> # Create a new branch and check it out
          The remote branch is automatically created when you push it to the remote server. So when you feel ready for it, you can just do:
  
          git push <remote-name> <branch-name>
          Where <remote-name> is typically origin, the name which git gives to the remote you cloned from. Your colleagues would then just pull that branch, and it's au    tomatically created locally.
  
          Note however that formally, the format is:
  
          git push <remote-name> <local-branch-name>:<remote-branch-name>
          But when you omit one, it assumes both branch names are the same. Having said this, as a word of caution, do not make the critical mistake of specifying only     :<remote-branch-name> (with the colon), or the remote branch will be deleted!
  
          So that a subsequent git pull will know what to do, you might instead want to use:
  
          git push --set-upstream <remote-name> <local-branch-name>
  
  
  Add a new file to the remote branch:
          git add .
          # Adds the file to your local repository and stages it for commit. To unstage a file, use 'git reset HEAD YOUR-FILE'.
          git commit -m "Add existing file"
          # Commits the tracked changes and prepares them to be pushed to a remote repository. To remove this commit and modify the file, use 'git reset --soft HEAD~1'     and commit and add the file again.
          git push origin your-branch
          # Pushes the changes in your local repository up to the remote repository you specified as the origin
