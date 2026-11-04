// Run with GroovyMain from the nf-test JAR; evaluate the guard used by the pipeline test.
def sourceFile = new File(args ? args[0] : 'tests/default.nf.test')
def source = sourceFile.getText('UTF-8')
def start = source.indexOf('def assertWorkflowVersions =')
def end = source.indexOf('def rawVersionsFile =', start)
assert start >= 0 && end > start, 'Cannot find the actual raw versions guard'
def guard = source.substring(start, end)
def check = new GroovyShell().evaluate(guard + '\nassertWorkflowVersions')

def positive = [
    release: "Workflow:\n  nf-core/atacseq: v2.2.0\n  Nextflow: 25.10.4\nFIXTURE_TOOL:\n  fixture: 1.2.3\n",
    development: "Workflow:\n  nf-core/atacseq: v2.2.0dev-gabcdef0\n  Nextflow: 26.10.5-edge\n",
]
def negative = [
    missing_workflow: 'FIXTURE_TOOL:\n  fixture: 1.2.3\n',
    missing_nextflow: 'Workflow:\n  nf-core/atacseq: v2.2.0dev\n',
    missing_pipeline_version: 'Workflow:\n  Nextflow: 25.10.4\n',
    other_pipeline: 'Workflow:\n  nf-core/other: v2.2.0\n  Nextflow: 25.10.4\n',
    duplicate_workflow: positive.release + 'Workflow:\n  nf-core/atacseq: v2.2.0\n  Nextflow: 25.10.4\n',
    duplicate_nextflow: 'Workflow:\n  nf-core/atacseq: v2.2.0\n  Nextflow: 25.10.4\n  Nextflow: 26.10.5\n',
    duplicate_pipeline_version: 'Workflow:\n  nf-core/atacseq: v2.2.0\n  nf-core/atacseq: v2.3.0\n  Nextflow: 25.10.4\n',
    malformed: 'Workflow: [',
    empty: '',
    empty_mapping: '{}\n',
    empty_workflow: 'Workflow: {}\n',
    nonmapping_workflow: 'Workflow: unexpected\n',
    nonmapping_document: '- Workflow\n',
    empty_nextflow: 'Workflow:\n  nf-core/atacseq: v2.2.0\n  Nextflow: ""\n',
    invalid_nextflow: 'Workflow:\n  nf-core/atacseq: v2.2.0\n  Nextflow: current\n',
    invalid_pipeline_version: 'Workflow:\n  nf-core/atacseq: unknown\n  Nextflow: 25.10.4\n',
]
def observations = []
positive.each { name, text ->
    check(text)
    observations << [case: name, expected: 'accept', observed: 'accept']
}
negative.each { name, text ->
    def rejected = false
    try {
        check(text)
    } catch (AssertionError | org.yaml.snakeyaml.error.YAMLException expected) {
        rejected = true
    }
    assert rejected, "Actual guard accepted negative control: ${name}"
    observations << [case: name, expected: 'reject', observed: 'reject']
}
def sha = { byte[] bytes -> java.security.MessageDigest.getInstance('SHA-256').digest(bytes).encodeHex().toString() }
def dependencyVersion = { String path ->
    def properties = new Properties()
    def stream = this.class.classLoader.getResourceAsStream(path)
    assert stream != null, "Missing runtime identity: ${path}"
    stream.withCloseable { properties.load(it) }
    properties.getProperty('version')
}
def receipt = [status: 'PASS', positive: positive.size(), negative: negative.size(), observations: observations,
    guard_sha256: sha(guard.getBytes('UTF-8')), test_source_sha256: sha(sourceFile.bytes),
    java_version: System.getProperty('java.version'),
    nf_test_version: dependencyVersion('META-INF/maven/com.askimed/nf-test/pom.properties'),
    snakeyaml_version: dependencyVersion('META-INF/maven/org.yaml/snakeyaml/pom.properties'),
    scope: 'Actual embedded guard evaluated with the nf-test runtime; version presence/syntax only, not a comparison to a separate actual-engine binding. No pipeline or snapshot execution.']
def json = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(receipt)) + '\n'
if (args.size() > 1) new File(args[1]).setText(json, 'UTF-8')
println json
