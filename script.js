const correctAnswers = {
    q1: 'b',
    q2: 'b',
    q3: 'c',
    q4: 'b',
    q5: 'a',
    q6: 'a',
    q7: 'c',
    q8: 'a',
    q9: 'a',
    q10: 'b'
  };
  
  function checkAnswers() {
    let score = 0;
    for (let i = 1; i <= 10; i++) {
      const selected = document.querySelector(`input[name="q${i}"]:checked`);
      const result = document.getElementById(`result${i}`);
      if (selected) {
        if (selected.value === correctAnswers[`q${i}`]) {
          result.textContent = "Correct!";
          result.className = "result correct";
          score++;
        } else {
          result.textContent = "Wrong!";
          result.className = "result incorrect";
        }
      } else {
        result.textContent = "Please select an answer.";
        result.className = "result incorrect";
      }
    }
    document.getElementById("finalScore").textContent = `Your total score is ${score}/10`;
  }
  