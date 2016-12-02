1. To add/delete file to github:
git add/rm  filename
2. To  make the commit
git commit -m 'I made those change becasue...'
3. Sync
git push
4. git clone https://....
5. git pull
6. git status


You might receive 403 error. There are many ways to solve it. One way is to input password for every push.

1. git config -l | grep url

Result: remote.origin.url=https://github.com/github-username/github-repository-name.git

2. git remote set-url origin "https://github-username@github.com/github-username/github-repository-name.git"

OR git remote set-url origin "https://Kee-Wang@github.com/Kee-Wang/repo.git"


## In order to avoid entering password everytime, one has to use ssh for each repo.
To do so:

1. Find the ssh rsh in currt node:

	` cd ~`

	`ssh-keygen -t rsa`    #Press enter for all values

2. Associate SSH to remote repo:

	Go to setting in the GitHub accout, and add ssh key, the key is in:

	`~/.ssh/id_rsa.pub`

3. Set a URL form that support SSH

	Desired form: `git+ssh://git@github.com/username/reponame.git`

	Wrong form: `https://github.com/username/reponame.git`

	In order to check: `git remote show origin`

	In oder to chang: `git remote set-url origin git+ssh://git@github.com/username/reponame.git`

	Or, in my case: 

	`git remote set-url origin git+ssh://git@github.com/Kee-Wang/repo.git`

4. You can restore a uncommit deletion by `git checkout -- `deleted-file-here`. Or you can follow the procedure by typing `git status`

5. How to cancel all current commits:

	git reset --hard origin/<branch_name>

	branch_name can be found: git branch
