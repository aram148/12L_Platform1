import os
import libcellml
#Parse a CellML file into a model m
p=libcellml.Parser()
with open('./experiment1/hai1.cellml') as f:
    contents=f.read()
    m=p.parseModel(contents)

print(m.componentCount())
print('The model has imports?',m.hasImports())
#print(m.componentCount())
# Create an Importer to resolve and flatten the model
i=libcellml.Importer()
result=i.resolveImports(m,'./experiment1/')
print('Resolving result:',result)

print('The importer found {} issues.'.format(i.issueCount()))
for index in range(i.issueCount()):
    print(i.issue(index).description())
	
model = i.flattenModel(m)

print(m.componentCount())
# Create a validator instance v and pass the model to it
v=libcellml.Validator()
# Validate the model m
v.validateModel(model)

a=libcellml.Analyser()
#Analyze the model
a.analyseModel(model)
	
print('The analyser found {} issues.'.format(a.issueCount()))
for b in range(0, a.issueCount()):
	issue = a.issue(b)
	print(issue.description())
		
# Retrieve the error information
print('The error count is',v.errorCount())
for index in range(v.errorCount()):
    print(v.error(index).description())

#print('The issue count',v.issueCount())
#for index in range(v.issueCount()):
#    print(v.issue(index).description())

#print('The warning count',v.warningCount())
#for index in range(v.warningCount()):
#    print(v.warning(index).description())

#print('The hint count',v.hintCount())
#for index in range(v.warningCount()):
#    print(v.hint(index).description())

#print('The message count',v.messageCount())
#for index in range(v.warningCount()):
#    print(v.message(index).description())

for d in range(0, i.libraryCount()):

        # Retrieve the library model by index, d.
    v.validateModel(i.library(d))

        # Retrieve the key under which it's stored: this will be the URL at which the imported model was found.
    print("The validator found {} issues in {}.".format(v.issueCount(),i.key(d)))
    for aa in range(0, v.issueCount()):
        print("    - {}".format(v.issue(aa).description()))

print()

f.close()